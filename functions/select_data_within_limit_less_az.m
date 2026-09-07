function [new,new_time]=select_data_within_limit_less_az(data,data_time,which_az)
    k_k=find((data(:,1)<which_az));
    new(:,:)=data(k_k,:); %Azimuth, Zenith, DoLP, AoP, DoCP
    if isempty(data_time)
    new_time=[];    
    else
    new_time(:,:)=data_time(k_k,:);
    end
end
