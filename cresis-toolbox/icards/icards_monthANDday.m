function [ month1,month2,day1 ] = icards_monthANDday( today )
today=num2str(today);
month_temp=today(5:6);

switch month_temp
  case '01'
    month1=1;month2='jan';
  case '02'
    month1=2;month2='feb';
  case '03'
    month1=3;month2='mar';
  case '04'
    month1=4;month2='april';
  case '05'
    month1=5;month2='may';
  case '06'
    month1=6;month2='jun';
  case '07'
    month1=7;month2='jul';  
  case '08'
    month1=8;month2='aug';
  case '09'
    month1=9;month2='sep';
  case '10'
    month1=10;month2='oct';  
  case '11'
    month1=11;month2='nov';  
  case '12'
    month1=12;month2='dec';
end

if str2num(today(7))==0
  day1=str2num(today(8));
else
  day1=str2num(today(7:8));
end


end

