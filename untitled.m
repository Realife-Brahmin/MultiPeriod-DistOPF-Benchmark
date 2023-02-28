% Input pairs of numbers as an array
pairs = [4 6; 2 5; 6 1; 4 5; 2 6; 1 5];

% Find the number of unique nodes in the pair array
nodes = unique(pairs);

% Create a mapping from original node number to new node number
mapping = zeros(max(nodes), 1);
mapping(nodes) = 1:length(nodes);

% Convert pairs to new node numbers
new_pairs = mapping(pairs);

% Create an adjacency matrix from the new pairs
adj_matrix = full(sparse(new_pairs(:, 1), new_pairs(:, 2), 1));

% Plot the graph
G = graph(adj_matrix);
plot(G);
