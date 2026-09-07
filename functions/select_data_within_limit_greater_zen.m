function [new,new_time]=select_data_within_limit_greater_zen(data,data_time,which_zenith)
    k_k=find((data(:,2)>which_zenith));
    new(:,:)=data(k_k,:); %Azimuth, Zenith, DoLP, AoP, DoCP
    if isempty(data_time)
    new_time=[];    
    else
    new_time(:,:)=data_time(k_k,:);
    end
end
