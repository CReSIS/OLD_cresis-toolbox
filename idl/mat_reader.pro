; function mat_reader (and support functions)
;
; Description:
; Load Matlab .MAT -v6 files
;
; EXAMPLE:
; From MATLAB:
; >> Data = [1 2 3 4]
; >> save('-v6','test.mat','Data')
; From IDL:
; .compile mat_reader.pro
; result = mat_reader('test.mat','Data', data)
;
; Author: John Paden

; ====================================================
; ====================================================
; mat_reader_data
;   Support function for mat_reader
; ====================================================
; ====================================================
function mat_reader_data, lun, tag, data, align=align

MAT_READER_miINT8 = 1
MAT_READER_miUINT8 = 2
MAT_READER_miINT16 = 3
MAT_READER_miUINT16 = 4
MAT_READER_miINT32 = 5
MAT_READER_miUINT32 = 6
MAT_READER_miSINGLE = 7
MAT_READER_miDOUBLE = 9
MAT_READER_miINT64 = 12
MAT_READER_miUINT64 = 13

CASE tag.data_type OF

MAT_READER_miINT8: data = bytarr(tag.num_bytes/1,/NOZERO)
MAT_READER_miUINT8: data = bytarr(tag.num_bytes/1,/NOZERO)
MAT_READER_miINT16: data = intarr(tag.num_bytes/2,/NOZERO)
MAT_READER_miUINT16: data = uintarr(tag.num_bytes/2,/NOZERO)
MAT_READER_miINT32: data = lonarr(tag.num_bytes/4,/NOZERO)
MAT_READER_miUINT32: data = ulonarr(tag.num_bytes/4,/NOZERO)
MAT_READER_miSINGLE: data = fltarr(tag.num_bytes/4,/NOZERO)
MAT_READER_miDOUBLE: data = dblarr(tag.num_bytes/8,/NOZERO)
MAT_READER_miINT64: data = lon64arr(tag.num_bytes/8,/NOZERO)
MAT_READER_miUINT64: data = ulon64arr(tag.num_bytes/8,/NOZERO)
ELSE: print,tag.data_type,format='(%"Unsupported type %d")'

ENDCASE

readu,lun,data

if NOT keyword_set(align) then begin
  align = 8
endif

; Make sure the file pointer aligns to the next (align)-byte boundary
point_lun,-lun,cur_pos
new_pos = cur_pos + ceil(double(tag.num_bytes)/align)*align - tag.num_bytes
if cur_pos NE new_pos then begin
  ;print,new_pos-cur_pos,format='(%"  Moving forward %d")'
  point_lun,lun,new_pos
endif

return,0

end

; ====================================================
; ====================================================
; mat_reader_tag
;   Support function for mat_reader
; ====================================================
; ====================================================
function mat_reader_tag, lun, tag, tag_data

tmp = lonarr(1)
readu,lun,tmp
if ishft(tmp,-16) EQ 0 then begin
  ;print, '  Regular data block'
  tag.data_type = tmp
  readu,lun,tmp
  tag.num_bytes = tmp
  tag.data_packed = 0
  data = [0]
endif else begin
  ;print, '  Data block with data in tag'
  tag.data_type = ishft(ishft(tmp,16),-16)
  tag.num_bytes = ishft(tmp,-16)
  tag.data_packed = 1
  result = mat_reader_data(lun, tag, tag_data, align=4)
endelse
;print, tag.data_type, format='(%"  data type: %d")'
;print, tag.num_bytes, format='(%"  num bytes: %d")'

return, 0

end

; ====================================================
; ====================================================
; mat_reader
;   file_name: string containing Matlab .MAT -v6 file name
;   var_name: string containing variable name that you
;     want to load (e.g. 'Data', 'Latitude', etc.)
;   data: Data matrix
;   Function returns 0 when variable is found
;     data is set to variable
;   Function returns 1 when variable is not found,
;     data is unchanged
; ====================================================
; ====================================================
function mat_reader, file_name, var_name, data

;print, file_name, format='(%"File name: %s")'

; ==============================================
; Open File
; ==============================================
openr,lun,file_name,/get_lun

; ==============================================
; Read header
; ==============================================
header = {text:bytarr(116), subsys_data_offset:bytarr(8), $
  version:intarr(1), endian_indicator:bytarr(2)}

readu,lun,header

;print, format='(%"Header")'
;print, string(header.text), format='(%" text: %s")'
;print, string(header.subsys_data_offset), format='(%" subsys_data_offset: %s")'
;print, string(header.version), format='(%" version: %d")'
;print, string(header.endian_indicator), format='(%" endian indicator: %s")'

; ==============================================
; ==============================================
; Read data blocks until end of file or data
; block of interest found
; ==============================================
; ==============================================
tag = {num_bytes:0L, data_type:0L, $
  data_packed:0}

var_num = 1
while (not EOF(lun)) do begin
  ; ==============================================
  ; Read in main variable tag
  ; ==============================================
  ;print, 'Variable' + string(var_num)
  result = mat_reader_tag(lun,tag,tag_data)

  point_lun,-lun,next_field
  next_field = next_field + tag.num_bytes

  if tag.data_type NE 14 then begin
    ;print, '  Variable must be miMATRIX type'
    point_lun,lun,next_field
    continue
  endif

  ; ==============================================
  ; Read in array flags subelement
  ; ==============================================
  ;print, ' Array flags'
  result = mat_reader_tag(lun,tag,tag_data)
  result = mat_reader_data(lun,tag,data)
  data_class = ishft(ishft(data[0],24),-24)
  ;print, data[0], format='(%"  flag 1: %d")'
  ;print, data[1], format='(%"  flag 2: %d")'
  ;print, data_class, format='(%"  data class: %d")'

  ; Check if this data class is supported
  if NOT (data_class EQ 4 OR (data_class GE 6 AND data_Class LE 15)) then begin
    ;print, '  Data class not supported'
    point_lun,lun,next_field
    continue
  endif

  ; ==============================================
  ; Read in dimensions subelement
  ; ==============================================
  ;print, ' Dimensions'
  result = mat_reader_tag(lun,tag,tag_data)
  result = mat_reader_data(lun,tag,dims)
  ;print, dims[0], dims[1], format='(%"  %d by %d")'

  ; ==============================================
  ; Read in array name subelement
  ; ==============================================
  ;print, ' Array Name'
  result = mat_reader_tag(lun,tag,tag_data)
  if tag.data_packed EQ 1 then begin
    array_name = tag_data;
  endif else begin
    result = mat_reader_data(lun,tag,array_name)
    array_name = string(array_name)
  endelse
  ;print, string(array_name), format='(%"  %s")'

  ; Check if this is the variable that we are loading
  if strcmp(var_name,string(array_name)) EQ 0 then begin
    ;print, var_name, string(array_name), format='(%"  Variable %s != %s")'
    point_lun,lun,next_field
    continue
  end

  ; ==============================================
  ; Read in data subelement
  ; ==============================================
  ;print, ' Data'
  result = mat_reader_tag(lun,tag,tag_data)
  if tag.data_packed EQ 1 then begin
    data = tag_data;
  endif else begin
    result = mat_reader_data(lun,tag,data)
  endelse

  data = reform(data,dims[0],dims[1])
  ;print, 'Found variable'
  free_lun,lun

  return,0

endwhile

;print, 'Did not find variable'

free_lun,lun

return, 1

end

