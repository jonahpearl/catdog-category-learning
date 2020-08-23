function [X, rIntervals] = load_interval_data(fname)
% load interval data and stitch if needed
Data = load(fname);
vars = fieldnames(Data);
if ismember('X1', vars)
    X = vertcat(Data.X1, Data.X2);
else
    X = Data.X;
end
rIntervals = Data.rIntervals;
