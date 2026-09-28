"""Read a .msh file, return {node_tag: (x, y, z)}."""
import sys
import gmsh


def get_node_coords(path):
    gmsh.initialize(["-v", "0"])
    gmsh.open(path)
    tags, xyz, _ = gmsh.model.mesh.getNodes()
    gmsh.finalize()
    xyz = xyz.reshape(-1, 3)
    return {int(t): tuple(xyz[i]) for i, t in enumerate(tags)}


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "./GFG.msh"
    nodes = get_node_coords(path)
    print(f"{len(nodes)} nodes")
    for tag in list(nodes)[:5]:
        print(tag, nodes[tag])
