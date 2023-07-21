% This script changes all interpretation parameters to LaTeX. Plot titles, labels, ticks, etc. will now be shown in the superior LaTeX representation instead of the generic, uninspiring regular MATLAB font.

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end