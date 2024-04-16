function total_length = total_length(data)
%% Calculates total length in um
    total_length = sum(data.Graph.segInfo.segLen_um(:));
end