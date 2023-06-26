function hdr_loc = get_hdr_loc(data)

index0 = find(data == hex2dec('AA'));
index1 = find(data(index0+1) == hex2dec('AA'));
index1 = index0(index1)+1;
index2 = find(data(index1+1) == hex2dec('55'));
index2 = index1(index2)+1;
index3 = find(data(index2+1) == hex2dec('55'));
index3 = index2(index3)+1;
index4 = find(data(index3+1) == hex2dec('AA'));
index4 = index3(index4)+1;
index5 = find(data(index4+1) == hex2dec('AA'));
index5 = index4(index5)+1;
%index6 = find(data(index5+1) == hex2dec('55'));
%index6 = index5(index6)+1;
%index7 = find(data(index6+1) == hex2dec('55'));
%index7 = index6(index7)+1;
hdr_loc = index5-5;
