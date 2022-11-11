"""
Automatically load tapes in the Qualstar Q40 and pull the files from the given list.

TAPES_FILE should point to a file produced by run_get_raw_files.m and map files
to tapes.

Author: Reece Mathews
"""
from subprocess import check_output
from collections import defaultdict, namedtuple
from contextlib import contextmanager
from pathlib import Path
import re
import os
import shutil

TAPES_FILE = 'tapes.txt'
OUR_TAPES = 114, 120, 122, 124, 208, 232
TAPE_LIB_DEV = "/dev/sg4"
TAPE_MOUNT_PATH = "/mnt/ltfs"

# Determine with `lsscsi -g`
TAPE_DRIVE_DEVS = {
    "0": "/dev/sg3",
    "1": "/dev/sg5",
}


def removeprefix(string, prefix):
    """Reimplementing the Python 3.9 method for removing a string prefix."""
    return string[len(prefix):] if string.startswith(prefix) else string

def removesuffix(string, suffix):
    """Reimplementing the Python 3.9 method for removing a string suffix."""
    return string[:-len(suffix)] if string.endswith(suffix) else string


def get_tape_num(tape):
    """Retrieve the number from the tape label."""
    return int(removesuffix(removeprefix(tape, "OIB"), "L8"))


def parse_tapes_file():
    """Read the tapes file into a dict."""
    tape_mapping = defaultdict(list)
    path_mapping = {}
    all_files = set()

    with open("tapes.txt") as f:
        season = None
        for line in f:
            if season is None:
                season = line
                continue
            if line.startswith("Filelist: "):
                continue
            if line.startswith("tapes filename stored_filename"):
                continue
            if line.strip() == "":
                season = None
                continue

            parts = line.split()
            tapes = parts[0].split(",")
            original_path = parts[1]
            tape_path = parts[2]

            path_mapping[tape_path] = original_path

            if tape_path in all_files:
                input("Duplicate file " + tape_path)
            all_files |= {tape_path}

            for tape in tapes:
                tape_num = get_tape_num(tape)

                if tape_num % 2 == 0:
                    tape_mapping[tape].append(tape_path)
                    break
            else:
                raise RuntimeError("No even tape for " + tape_path)

    return tape_mapping, path_mapping


def inventory():
    """Check the tape library's current inventory."""
    INV_RE_MATCH_INDICES = {
        "slot_type": 0,
        "slot_num": 1,
        "mailslot": 2,
        "full": 3,
        "original_slot_num": 6,
        "barcode": 7
    }
    # Hope someone else doesn't end up working on this 👍  - reece
    INV_RE_PATTERN = r"(Storage|Data Transfer) Element ([0-9]+)\s?(IMPORT\/EXPORT)?:(Full|Empty)( (\(Storage Element ([0-9]+) Loaded\))?:VolumeTag\s?=\s?(\w+))?"
    InvSlot = namedtuple("InvSlot", "slot_type slot_num mailslot full original_slot_num barcode")

    inv_output = check_output(["mtx", "-f", TAPE_LIB_DEV, "status"]).decode()
    inv_output = "\n".join(line.strip() for line in inv_output.split("\n"))
    matches = re.findall(INV_RE_PATTERN, inv_output)
    return [InvSlot(**{k: match[i] for k, i in INV_RE_MATCH_INDICES.items()}) for match in matches]


def load_tape(tape_barcode):
    """Load a tape into a drive"""
    inv = inventory()

    # Find tape
    slot = None
    for slot_ in inv:
        if slot_.barcode == tape_barcode:
            slot = slot_
            break
    else:
        raise RuntimeError("Tape not present: " + tape_barcode)

    if slot.slot_type == "Data Transfer":
        # tape already in drive
        return slot.slot_num

    # Find open drive
    drive_num = None
    for slot_ in inv:
        if slot_.slot_type == "Data Transfer":
            drive_num = slot_.slot_num
            if slot_.full == "Empty":
                break
    else:
        # All drives full, unload last drive found
        if drive_num is None:
            raise RuntimeError("No drives found")

        check_output(["mtx", "-f", TAPE_LIB_DEV, "unload", slot.original_slot_num, str(drive_num)])

    check_output(["mtx", "-f", TAPE_LIB_DEV, "load", slot.slot_num, drive_num])

    return drive_num


def unmount_drive():
    """Unmount the TAPE_MOUNT_PATH."""
    if os.path.ismount(TAPE_MOUNT_PATH):
        check_output(["umount", TAPE_MOUNT_PATH])


@contextmanager
def mount_drive(drive_num):
    """Mount the given drive with LTFS."""
    unmount_drive()
    drive_dev = TAPE_DRIVE_DEVS[drive_num]
    check_output(["ltfs", "-o", f"devname={drive_dev}", TAPE_MOUNT_PATH])
    try:
        yield
    finally:
        unmount_drive()


def copy_files(files, path_mapping):
    """Copy all the given files to their original locations."""
    for file in files:
        # Perform sanity checks on file and destination
        # fs = filesystem (destination) paths as opposed to tape (source) paths

        file_path = TAPE_MOUNT_PATH + file
        file_exists = os.path.exists(file_path)
        file_size = os.path.getsize(file_path) / 1024 ** 2 if file_exists else None
        print(file_exists, f"{file_size}MB", file_path)

        fs_file_path = path_mapping[file]
        fs_file_exists = os.path.exists(fs_file_path)
        fs_file_size = os.path.getsize(fs_file_path) / 1024 ** 2 if fs_file_exists else None
        print("->", fs_file_exists, f"{fs_file_size}MB", fs_file_path)

        # Check if source exists on tape
        if not file_exists:
            input("**Could not find file on tape**")
            continue

        # Check if destination already exists on filesystem
        if fs_file_exists:
            if input("**File already exists on file system**, skip (y) or halt (n)?") == "y":
                continue
            else:
                raise RuntimeError("File already exists on file system")

        # Check if parent folder path exists
        fs_parent_path = Path(fs_file_path).parent
        if not fs_parent_path.exists():
            if input("**Path to file does not exist**, create path (y) or skip (n)?") == "y":
                os.makedirs(fs_parent_path, exist_ok=True)
            else:
                continue

        # Perform copy
        shutil.copy2(file_path, fs_file_path)


def load_tapes(tape_mapping, path_mapping):
    """Load each tape and retrieve the files back to their original location."""
    for tape in tape_mapping:

        try:
            drive_num = load_tape(tape)
        except RuntimeError:
            print("Skipping missing tape: " + tape)
            continue

        files = tape_mapping[tape]
        with mount_drive(drive_num):
            copy_files(files, path_mapping)


if __name__ == "__main__":
    if os.geteuid() != 0:
        raise RuntimeError("Must be ran as root to access tape device")

    tape_mapping, path_mapping = parse_tapes_file()
    load_tapes(tape_mapping, path_mapping)
