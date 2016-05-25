function result = strip_hdr(fin,fout)

%rewind the files
fseek(fin,0,-1);
fseek(fout,0,-1);

%find first file sync
data = fread(fin,10000,'uint8');
hdr_loc = hfrds.get_hdr_loc(data);

%read header
[hdr,check] = hfrds.get_hdr(fin,hdr_loc(1)-1);

%continue reading new data
i=0;
while(check),
  %write hdr data to a file
  hfrds.write_hdr(fout,hdr);
  i=i+1; epri = hdr.epri;
  %read next header
  [hdr,check] = hfrds.next_hdr(fin,hdr);
  %make sure epri is counting by ones (report error if not)
  if (hdr.epri-epri) ~= 1, check=0; end;
end;

result = i;
