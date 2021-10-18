import bisect
import copy
import timeit

class DeBruijnGraph:
    def __init__(self, text: str, k: int = 3):
        self.text = text
        self.k = k
        self.graph_edges = []
        self.graph_nodes = {}  # values are lists of 2 ints: [number of
        # occurences on left side of k-mer, number of occurences on right side of k-mer]
        self.prepare()

    def prepare(self):
        for i in range(len(self.text) - self.k + 1):
            self.graph_edges.append((self.text[i:i + self.k - 1], self.text[i + 1:i + self.k]))
            self.graph_nodes.setdefault(self.text[i:i + self.k - 1], [0, 0])[0] += 1
            self.graph_nodes.setdefault(self.text[i + 1:i + self.k], [0, 0])[1] += 1


def find_first_node(nodes):
    try:
        first_node = [node for node, occurences in nodes.items()
                      if occurences[1] + 1 == occurences[0]][0]
        return first_node
    except:
        return None

def find_last_node(nodes):
    last_node = [node for node, occurences in nodes.items()
                 if occurences[0] + 1 == occurences[1]][0]
    return last_node

def left_right_index(edges, node):
    index_left = bisect.bisect_left(edges, node)
    index_right = bisect.bisect_right(edges, node)
    return index_left, index_right

def eulerian_walk_with_two_bisects(nodes, edges, reconstructed_text, list_of_walks):
    current_node = find_first_node(nodes)
    if current_node == None:
        if len(edges) == 0:
            list_of_walks.add(reconstructed_text)
        return list_of_walks

    last_node = find_last_node(nodes)
    while True:
        nodes[current_node][0] = nodes[current_node][0] - 1

        if current_node == last_node and nodes[last_node][0] == -1:
            if len(edges) == 0:
                list_of_walks.add(reconstructed_text)
            return list_of_walks

        index_left = bisect.bisect_left([edge[0] for edge in edges], current_node)
        index_right = bisect.bisect_right([edge[0] for edge in edges], current_node)

        if index_right - index_left > 1:
            possible_edges_set = set(edges[index_left:index_right])
            for i in range(index_left, index_right):
                if edges[i] in possible_edges_set:
                    current_node = edges[i][1]
                    eulerian_walk_kwargs = {'nodes': copy.deepcopy({key: (value if key != current_node
                                                                    else [value[0], value[1] - 1])
                                                                    for key, value in nodes.items()}),
                                            'edges': copy.deepcopy(edges[:i] + edges[i+1:]),
                                            'reconstructed_text': reconstructed_text + current_node[-1],
                                            'list_of_walks': list_of_walks}
                    list_of_walks = eulerian_walk_with_two_bisects(**eulerian_walk_kwargs)
                    possible_edges_set.remove(edges[i])
            return list_of_walks
        else:
            current_node = edges.pop(index_left)[1]
            nodes[current_node][1] -= 1
            reconstructed_text += current_node[-1]

def eulerian_walk_with_index_count_func(nodes, edges, reconstructed_text, list_of_walks):
    current_node = find_first_node(nodes)
    if current_node == None:
        if len(edges) == 0:
            list_of_walks.add(reconstructed_text)
        return list_of_walks

    last_node = find_last_node(nodes)
    while True:
        nodes[current_node][0] = nodes[current_node][0] - 1

        if current_node == last_node and nodes[last_node][0] == -1:
            if len(edges) == 0:
                list_of_walks.add(reconstructed_text)
            return list_of_walks

        index_left, index_right = left_right_index([edge[0] for edge in edges], current_node)

        if index_right - index_left > 1:
            possible_edges_set = set(edges[index_left:index_right])
            for i in range(index_left, index_right):
                if edges[i] in possible_edges_set:
                    current_node = edges[i][1]
                    eulerian_walk_kwargs = {'nodes': copy.deepcopy({key: (value if key != current_node
                                                                    else [value[0], value[1] - 1])
                                                                    for key, value in nodes.items()}),
                                            'edges': copy.deepcopy(edges[:i] + edges[i+1:]),
                                            'reconstructed_text': reconstructed_text + current_node[-1],
                                            'list_of_walks': list_of_walks}
                    list_of_walks = eulerian_walk_with_index_count_func(**eulerian_walk_kwargs)
                    possible_edges_set.remove(edges[i])
            return list_of_walks
        else:
            current_node = edges.pop(index_left)[1]
            nodes[current_node][1] -= 1
            reconstructed_text += current_node[-1]

if __name__ == "__main__":
    setup = '''
    DB_graph = DeBruijnGraph("to_every_thing_turn_turn_turn_there_is_a_season", 3)
    edges = copy.deepcopy(DB_graph.graph_edges)
    nodes = copy.deepcopy(DB_graph.graph_nodes)
    eulerian_walk_kwargs = {'nodes': dict(sorted(nodes.items())),
                            'edges': sorted(edges),
                            'reconstructed_text': find_first_node(nodes),
                            'list_of_walks': set()}
    '''
    code_to_execute = '''eulerian_walk_with_index_count_func(**eulerian_walk_kwargs)'''
    # for i, walk in enumerate(list_of_walks):
    #     print(i, walk)
    one_list_comprehension = timeit.repeat(stmt=code_to_execute.replace('    ', ''),
                                           setup=setup.replace('    ', ''), globals=globals(),
                                           repeat=30)
    print('One list comprehension: ', one_list_comprehension, ' Average: ',
          sum(one_list_comprehension)/len(one_list_comprehension))

    code_to_execute = '''eulerian_walk_with_two_bisects(**eulerian_walk_kwargs)'''
    two_list_comprehension = timeit.repeat(stmt=code_to_execute.replace('    ', ''),
                                           setup=setup.replace('    ', ''), globals=globals(),
                                           repeat=30)
    print('Two list comprehension: ', two_list_comprehension, ' Average: ',
          sum(two_list_comprehension) / len(two_list_comprehension))



