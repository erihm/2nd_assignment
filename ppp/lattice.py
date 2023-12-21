import itertools
import graphviz

def power_set(s):
    power_set_list = []
    for i in range(len(s) + 1):
        power_set_list.extend(itertools.combinations(s, i))
    return [set(subset) for subset in power_set_list]

def draw_power_set_lattice(subsets_to_plot):
    graph = graphviz.Digraph('PowerSetLattice', format='png', graph_attr={'rankdir': 'BT'})

    for i, subset in enumerate(subsets_to_plot):
        graph.node(str(i), ', '.join(map(str, subset)))

    for i, subset in enumerate(subsets_to_plot):
        for j, other_subset in enumerate(subsets_to_plot):
            if subset.issubset(other_subset) and i != j and len(subset) + 1 == len(other_subset):  # Check for immediate superset
                graph.edge(str(i), str(j))

    graph.render(filename='power_set_lattice', format='png', cleanup=True)

if __name__ == "__main__":
    set_elements = {1, 2, 3} 
    subsets_to_plot = [{2}, {3}, {1, 2}, {2, 3}, {1, 2, 3}]
    draw_power_set_lattice(subsets_to_plot)
