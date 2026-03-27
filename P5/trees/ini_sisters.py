from tree_nodes import Node, newick2nodes


def sisters(nodes):
    """
    Returns a list of sister leaves. Each sister leaf is represented as a tuple: left child, right child, and parent node ids.

    >>> sisters(newick2nodes("((A,B),((C,(D,E)),F));"))
    [(2, 3, 1), (8, 9, 7)]
    >>> sisters(newick2nodes("((A,B),(((C,D),E),F));"))
    [(2, 3, 1), (7, 8, 6)]
    >>> sisters(newick2nodes("(((hmgb_chite:0.10,hmgl_wheat:0.25):0.20,hmgl_trybr:0.60):0.25,hmgt_mouse:0.35);"))
    [(3, 4, 2)]
    """
    
    all_leaves = [
        (i, value) 
        for i, value in nodes.items()
        if value.right == 0
    ]
    all_sisters = []
    
    for i in range(len(all_leaves) - 1):
        if all_leaves[i][1].parent == all_leaves[i+1][1].parent:
            all_sisters.append((
                all_leaves[i][0],
                all_leaves[i+1][0],
                all_leaves[i][1].parent
                )
            )

    return all_sisters

if __name__ == "__main__":
    print(sisters(newick2nodes("((A,B),((C,(D,E)),F));")))
    print(
        sisters(
            newick2nodes(
                "(((hmgb_chite:0.10,hmgl_wheat:0.25):0.20,hmgl_trybr:0.60):0.25,hmgt_mouse:0.35);"
            )
        )
    )
