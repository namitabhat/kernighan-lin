% implementation of kernighan lin algorithm

% load edgelist from csv file and create adjacency matrix
% each line contains edge from node in the 1st column to node in the 2nd column
E = load('path to data file'); 
nr_edges = length(E);
N = max(max(E));
A = sparse(E(:,1),E(:,2),1,N,N);


% start with first half of nodes in group1 and second half of nodes in group2 
N1 = N/2;
N2 = N/2;
g1 = 1:N1;
g2 = N1+1:N;


while 1
    i = 1;
    
    % initialise forbidden set C
    C = [];
    
    % calculate nominal cut size Rab_nom
    g1nom = g1;
    g2nom = g2;
    Rabnom = calculate_Rab(A,g1nom,g2nom);
    
    % initialise Rab_swap, storing change in cut size after each iteration 
    % initialise g1_state and g2_state, storing elements of group1 and group2 after each iteration
    Rab_swap = zeros(1,N1);
    g1_state = zeros(N1,N1);
    g2_state = zeros(N1,N1);
    
    while (i<=N1)
        % find node pair (a,b) with largest reduction in cut size
        % a is a node in group1, b is a node in group 2
        max_diff = -99;
        for a=1:length(g1)
            for b=1:length(g2)
                g1([1 a]) = g1([a 1]);
                g2([1 b]) = g2([b 1]);
        
                % check a and b are not in forbidden set 
                if (~ismember(g1(1),C) && ~ismember(g2(1),C) && reduction(A,g1,g2)>max_diff)
                    a_largest = g1(1);
                    b_largest = g2(1);
                    max_diff = reduction(A,g1,g2);
                end
                g1([a 1]) = g1([1 a]);
                g2([b 1]) = g2([1 b]);
            end
        end
        
        % swap node pair (a_largest,b_largest) and add them to C
        % calculate cut size after swapping nodes 
        g1(g1==a_largest) = b_largest;
        g2(g2==b_largest) = a_largest;
        Rab_swap(1,i) = calculate_Rab(A,g1,g2);
        C = [C, a_largest, b_largest];
        
        %save current state of group1 and group2
        g1_state(i,:) = g1;
        g2_state(i,:) = g2;
        
        i = i + 1; 
    end
    
    % find iteration j with smallest cut size after swap
    j = find(Rab_swap==min(Rab_swap),1);
    if Rab_swap(1,j)>=Rabnom
        break;
    end
    
    % restart loop with the state in iteration j
    g1 = g1_state(j,:);
    g2 = g2_state(j,:);   
end

print_graph(A,g2nom,'kl_partition.gv');


% calculate reduction in cut size
function diff = reduction(A,g1,g2)

    a = g1(1);
    b = g2(1);
    
    a_row = A(a,:);
    b_row = A(b,:);
    Ea = sum(a_row(g2));
    Eb = sum(b_row(g1));
    Ia = sum(a_row(g1));
    Ib = sum(b_row(g2));
    
    diff = Ea + Eb - Ia - Ib - 2*A(a,b);
    
end


% calculate cut size
function Rab = calculate_Rab(A,g1,g2)

    a = g1(1);
    b = g2(1);
    
    a_row = A(a,:);
    b_row = A(b,:);
    Ea = sum(a_row(g2));
    Eb = sum(b_row(g1));
    
    Z = 0;
    g1(1) = [];
    g2(1) = [];
    for i=1:length(g1)
        for j=1:length(g2)
            Z=Z+A(g1(i),g2(j));
        end
    end
    Rab = Z+Ea+Eb-A(a,b);
end


% generate .gv file to plot partition graph 
function print_graph(adjacency_matrix,g2,filename)
    
% node_states specifies partition
% [0, 0, 1] means nodes 1 and 2 belong to group1 and node 3 belongs to group2
nr_nodes=size(adjacency_matrix,1);
node_states = zeros(1,nr_nodes);
node_states(:,g2) = 1;
labels ={'A','B','C','D'};
colors={'red','blue','black','green'};

file1=fopen(filename,'w');
fprintf(file1,'graph G {\n');

% define nodes
for current_node=1:nr_nodes
    current_node_label = ['"' num2str(current_node) '\n' labels{node_states(current_node)+1}  '"'  ' [color = ' colors{node_states(current_node)+1} ']'];
    fprintf(file1,'%s ; \n',current_node_label);
end

% connect nodes
for to_node=1:nr_nodes
    for from_node=to_node+1:nr_nodes
        if adjacency_matrix(to_node,from_node)
            label1 = ['"' num2str(to_node) '\n' labels{node_states(to_node)+1} '"'];
            label2 = ['"' num2str(from_node) '\n' labels{node_states(from_node)+1} '"'];
            fprintf(file1,'%s -- %s ; \n',label1,label2);
        end
    end
end
fprintf(file1,'}');
fclose(file1);

end



