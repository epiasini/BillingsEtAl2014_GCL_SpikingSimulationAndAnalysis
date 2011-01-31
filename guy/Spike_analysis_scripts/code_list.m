% Function to find hexadecimal code list from
% binary spike train data
% Guy Billings, UCL 2010

function list=code_list(chunked_data,smode,slicedex)

data_dims=size(chunked_data);
codes={};
for observation=1:data_dims(1)
        
        data=squeeze(chunked_data(observation,:,:,:));
        this_vec=slice(data,smode);
        hex=binar2hex(this_vec(slicedex,:));
        codes=[codes hex];
        
end
list=union(codes,codes);

% Convenient to assume that every code is unique, warn if not

ic=size(codes);
ic=ic(2);
uc=size(list);
uc=uc(2);

if ic~=uc
    warning(['Repeat input codes found, equiprobability assumption violated: ',num2str(ic-uc),' repeats'])
end    
