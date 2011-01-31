% Function to convert a binary array into a hexadecimal string 

function hex=binar2hex(binar)

len=size(binar);
len=len(2);
hex=[];

loop=floor(len/4);

for l=1:loop
    
    nstr=[num2str(binar(4*(l-1)+1)),num2str(binar(4*(l-1)+2)),num2str(binar(4*(l-1)+3)),num2str(binar(4*(l-1)+4))];
    num=bin2dec(nstr);
    hex=[hex dec2hex(num)];
    
end

read=loop*4;

if read~=len
    nstr=[];
    for l=1:len-read
        nstr=[nstr num2str(binar(read+l))];
    end    
    num=bin2dec(nstr);
    hex=[hex dec2hex(num)];
    
end    