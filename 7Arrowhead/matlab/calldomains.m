function [wut] = calldomains(binMatrFname, filteredMatrFname)
% Function to find connected components in a binary matrix, and then find the
% index of the highest value (local maximum) within that given connected
% component

% read data from file
binMatr = dlmread(binMatrFname);
filteredMatr = dlmread(filteredMatrFname);

% get connected components
[L, num] = bwlabel(binMatr);

domain_boundaries_r = []
domain_boundaries_c = []

% iterate through every connected component
for ccnum = 1:1:num
    % get the indices of pixels in that component
    [r, c] = find(L=ccnum);

    % for every idex
        % check the corresponding value in the the filtered matrix
        % if the value is larger than max, set max = index, note the index
    for index = 1:1:numel(r)
        tmp = filteredMatr(r(index), c(index));
        if index == 1 || tmp > maxVal
            maxVal = tmp;
            maxR = r(index);
            maxC = c(index);
        end
    end

    % save the index of the max value
    domain_boundaries_r(end + 1) = maxR;
    domain_boundaries_c(end + 1) = maxC;
end

% write the list of indices to a file, or return them as output
wut = [domain_boundaries_r domain_boundaries_c]
