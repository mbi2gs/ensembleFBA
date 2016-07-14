function [p1,p2] = uniquePerms(N)

p1 = randperm(N);
p2 = randperm(N);

% Make sure p1 and p2 are different
while sum(p1 == p2) == length(p1)
    p2 = randperm(N);
end

end