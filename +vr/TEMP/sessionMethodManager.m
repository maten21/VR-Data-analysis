
clear

path1 = 'C:\Users\Mate Neubrandt\Documents\MATLAB\NANSEN\code\integrations\sessionmethods';
path2 = 'C:\Users\Mate Neubrandt\UIO Physiology Dropbox Dropbox\Lab Data\Mate Neubrandt\NANSEN-project\mate_vr\Session Methods\+mate_vr';

d1 = dir([path1, '\**\*']);
dd1 = d1(~[d1.isdir]);

d2 = dir([path2, '\**\*']);


dd = vertcat(dd1,dd2);

dd2 = d2(~[d2.isdir]);
pathLength = length(path2);
for i = 1:1:height(dd2)
    dd2(i).folder = dd2(i).folder(pathLength+1:end);
end
table(3,4)
t = struct2table(dd2);
t = t(:,[1 2 5]);

t.Properties.VariableNames = {'Name','Menu location','Installed'};
t2 = table('Size',[height(t),2],'VariableTypes',{'string','string'});
t2.Properties.VariableNames = {'Description','Author'};

tt = [t t2];

fig = uifigure("Position",[20 20 650 650]);
uit = uitable(fig,"Data",tt,"Position",[20 20 600 600]);
uit.ColumnEditable = [false false true];
uit.ColumnSortable = true;

uit.Data.isdir




