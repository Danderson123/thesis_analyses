import networkx as nx
import json
import sys

def load_gml_graph(gml_file):
    """Load a GML file, make it a simple undirected graph with no duplicate edges."""
    G = nx.read_gml(gml_file)

    # Convert MultiDiGraph or DiGraph → undirected simple Graph
    if isinstance(G, (nx.MultiDiGraph, nx.MultiGraph)):
        undirected = nx.Graph()
        for u, v in G.edges():
            undirected.add_edge(u, v)
        return undirected
    else:
        return G.to_undirected()

def write_graph_to_gml(G, output_path):
    """Save the graph to a GML file for debugging."""
    nx.write_gml(G, output_path)
    print(f"Graph written to {output_path} for inspection.")

def extract_unitigs(G):
    """
    Extract unitigs (maximal non-branching paths) from an undirected graph.
    A unitig is a path where all internal nodes have degree 2.
    """
    visited = set()
    unitigs = []

    for node in G.nodes():
        if node in visited or G.degree(node) != 1:
            continue  # Start only from tips

        path = [node]
        visited.add(node)
        current = node
        prev = None

        while True:
            neighbors = [n for n in G.neighbors(current) if n != prev]
            if len(neighbors) != 1:
                break  # Dead end or branch point
            next_node = neighbors[0]
            if G.degree(next_node) != 2:
                path.append(next_node)
                visited.add(next_node)
                break  # End of unitig
            path.append(next_node)
            visited.add(next_node)
            prev = current
            current = next_node

        if len(path) > 1:
            unitigs.append(path)

    # Now look for cycles (all degree 2 and unvisited)
    for node in G.nodes():
        if node in visited:
            continue
        if G.degree(node) == 2:
            cycle = []
            current = node
            prev = None
            while True:
                cycle.append(current)
                visited.add(current)
                neighbors = [n for n in G.neighbors(current) if n != prev]
                if not neighbors:
                    break
                next_node = neighbors[0]
                if next_node == node:
                    cycle.append(node)  # complete cycle
                    break
                if next_node in visited:
                    break
                prev = current
                current = next_node
            if len(cycle) > 1:
                unitigs.append(cycle)

    return unitigs

def save_unitigs_to_json(unitigs, output_file):
    with open(output_file, 'w') as f:
        json.dump(unitigs, f, indent=2)

def main(gml_file, unitigs_json, debug_gml_out):
    G = load_gml_graph(gml_file)
    write_graph_to_gml(G, debug_gml_out)  # Save for debugging
    unitigs = extract_unitigs(G)
    save_unitigs_to_json(unitigs, unitigs_json)
    print(f"Extracted {len(unitigs)} unitigs and saved to {unitigs_json}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python extract_unitigs_from_gml.py <input.gml> <unitigs.json> <debug_out.gml>")
        sys.exit(1)

    gml_path = sys.argv[1]
    json_path = sys.argv[2]
    debug_gml_path = sys.argv[3]
    main(gml_path, json_path, debug_gml_path)
