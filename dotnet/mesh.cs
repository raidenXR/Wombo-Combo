namespace GridGeneration;

public struct Vec2
{
	public double X;
	public double Y;

	public Vec2(double x, double y)	{
		X = x;
		Y = y;
	}
}

public struct BoundingBox
{
	public Vec2 v_min;
	public Vec2 v_max;

}

public struct Index
{
	public int I;
	public int J;
}

public struct Gradient
{
	public double UL;
	public double UC;
	public double UR;
	public double CL;
	public double CR;
	public double LL;
	public double LC;
	public double LR;

	public unsafe GradDirection SteepestDescend()
	{
		Gradient copy = this;
		double* ptr  = (double*)&copy;
		int idx = 0;

		for (int i = 0; i < 8; i++) {
			idx = ptr[idx] > ptr[i] ? i : idx;
		}

		switch (idx) {
			case 0: return GradDirection.UL;
			case 1: return GradDirection.UC;
			case 2: return GradDirection.UR;
			case 3: return GradDirection.CL;
			case 4: return GradDirection.CR;
			case 5: return GradDirection.LL;
			case 6: return GradDirection.LC;
			case 7: return GradDirection.LR;
			default: return GradDirection.Self;
			
		}
	}

	
	public unsafe GradDirection SteepestAscend()
	{
		Gradient copy = this;
		double* ptr  = (double*)&copy;
		int idx = 0;

		for (int i = 0; i < 8; i++) {
			idx = ptr[idx] < ptr[i] ? i : idx;
		}

		switch (idx) {
			case 0: return GradDirection.UL;
			case 1: return GradDirection.UC;
			case 2: return GradDirection.UR;
			case 3: return GradDirection.CL;
			case 4: return GradDirection.CR;
			case 5: return GradDirection.LL;
			case 6: return GradDirection.LC;
			case 7: return GradDirection.LR;
			default: return GradDirection.Self;
			
		}
	}
}

public enum GradDirection
{
	UL, UC, UR, CL, CR, LL, LC, LR, Self,
}

public enum NodeState
{
	Internal, External,	Unknown,
}

public enum Marching
{
	Vertical, Horizontal,
}

// it holds only functions, no fields
public static class Geo
{
	public static Vec2 Center (Vec2 a, Vec2 b)
	{
		var x = a.X + (b.X - a.X) / 2.0;
		var y = a.Y + (b.Y - a.Y) / 2.0;

		return new Vec2(x,y);
	}

	public static Vec2[] Centers (Vec2[] vertices)
	{
		Vec2[] centers = new Vec2[vertices.Length-1];

		for (int i = 0; i < vertices.Length-1; i++) {
			var a = vertices[i+0];
			var b = vertices[i+1];
			var c = Center(a,b);
			centers[i] = c;
		}

		return centers;
	}

	public static Vec2 Normal (Vec2 a, Vec2 b)
	{
		var dx = b.X - a.X;
		var dy = b.Y - a.Y;
		var l  = Math.Sqrt(dx*dx + dy*dy);

		return new Vec2(dy, -dx);
	}

	public static Vec2[] Normals (Vec2[] vertices)
	{
		Vec2[] normals = new Vec2[vertices.Length-1];

		for (int i = 0; i < vertices.Length-1; i++) {
			var a = vertices[i+0];
			var b = vertices[i+1];
			var c = Center(a,b);
			var n = Normal(a,b);
			n.X += c.X;
			n.Y += c.Y;
			normals[i] = n;
		}

		return normals;
	}
	
	public static Vec2 Tangent (Vec2 a, Vec2 b)
	{
		return new Vec2(b.X - a.X, b.Y - a.Y);
	}

	/// <summary> transforms to x,y (Cartensian system) from i,j (Stencil system) </summary> 
	public static Vec2 ToCartesianSystem (bool[] grid, int N, int i, int j, BoundingBox bounds)
	{
		double x_min = bounds.v_min.X;
		double y_min = bounds.v_min.Y;
		double dx = (bounds.v_max.X - bounds.v_min.X) / (double)N;
		double dy = (bounds.v_max.Y - bounds.v_min.Y) / (double)N;		

		return new Vec2{
			X = (double)j * dx + x_min,
			Y = (double)i * dy + y_min,
		};		
	}

	/// <summary> transforms to i,j (Stencilsystem) from x,y (Cartesian system) </summary> 
	public static Index ToStencilSystem (int N, Vec2 p, BoundingBox bounds)
	{
		double x_min = bounds.v_min.X;
		double y_min = bounds.v_min.Y;
		double dx = (bounds.v_max.X - x_min) / (double)N;
		double dy = (bounds.v_max.Y - y_min) / (double)N;		
		
		int i = (int)Math.Round((p.Y - y_min) / dy, 0);
		int j = (int)Math.Round((p.X - x_min) / dx, 0);

		return new Index{
			I = Math.Clamp(i, 0, N-1),
			J = Math.Clamp(j, 0, N-1),
		};		
	}

	/// <summary> -- RECURSIVE FUNCTION --
	/// stencil method-algorithm for setting to 1 (true) the elements that are part of the boundary
	/// of the domain </summary>
	public static void AssignStencilElement (bool[] grid, int N, Vec2 a, Vec2 b, BoundingBox bounds)
	{
		var index_a = ToStencilSystem(N, a, bounds);
		var index_b = ToStencilSystem(N, b, bounds);

		// subdived (recursively) if nessacary
		if (Math.Abs(index_a.I - index_b.I) > 1 || Math.Abs(index_a.J - index_b.J) > 1) {
			AssignStencilElement(grid, N, a, Center(a,b), bounds);
			AssignStencilElement(grid, N, Center(a,b), b, bounds);
		}		

		grid[index_a.I * N + index_a.J] = true;
		grid[index_b.I * N + index_b.J] = true;
	}

	public static BoundingBox GetBounds (Vec2[] vertices)
	{
		int len = vertices.Length;
		double x_min = vertices[0].X;
		double x_max = vertices[0].X;
		double y_min = vertices[0].Y;
		double y_max = vertices[0].Y;
		
		for (int i = 1; i < len; i++) {
			x_min = x_min > vertices[i].X ? vertices[i].X : x_min;
			x_max = x_max < vertices[i].X ? vertices[i].X : x_max;
			y_min = y_min > vertices[i].Y ? vertices[i].Y : y_min;
			y_max = y_max < vertices[i].Y ? vertices[i].Y : y_max;
		}

		return new BoundingBox{
			v_min = new Vec2(x_min, y_min),
			v_max = new Vec2(x_max, y_max),		
		};
	}

	public static bool[] StencilFromBoundaryNodes (Vec2[] boundary_nodes, int N)
	{
		int len = boundary_nodes.Length;
		var bounds = GetBounds(boundary_nodes);
		var grid = new bool[N * N];
		
		for (int i = 0; i < len-1; i++) {
			Vec2 a = boundary_nodes[i+0];
			Vec2 b = boundary_nodes[i+1];
			AssignStencilElement(grid, N, a, b, bounds);
			
			Vec2 c = Center(a,b);
			Index ci = ToStencilSystem(N, c, bounds);
			grid[ci.I*N + ci.J] = true;
		}		

		return grid;
	}

	public static void FillStencil (bool[] grid, int N, BoundingBox bounds, Marching marching)
	{
		// for (int i = 0; i < N; i++) {
		// 	bool fill = grid[i*N + 0];

		// 	for (int j = 1; j < N-1; ) {
		// 		if (grid[i*N + j] && !fill) {
		// 			while (grid[i*N + j]) j++;
		// 			fill = true;
		// 		}				
		// 		else if (grid[i*N + j] && fill) {
		// 			while (grid[i*N + j]) j++;
		// 			fill = false;
		// 		}				
		// 		else if (fill) {
		// 			grid[i*N + j] = true;
		// 			j++;
		// 		}
		// 		else {
		// 			j++;
		// 		}
		// 	}
		// }
		switch (marching) {
			case Marching.Horizontal:
				for (int i = 0; i < N; i++) {
					int lhs = 0;
					int rhs = N-1;

					while (!grid[i*N + lhs]) lhs++; 
					while (!grid[i*N + rhs]) rhs--;

					while (grid[i*N + lhs]) lhs++;
					while (grid[i*N + rhs]) rhs--;

					int a = Math.Min(lhs, rhs);
					int b = Math.Max(lhs, rhs);
					bool fill = true;
					
					for (int j = a; j <= b; j++) {
						// if (grid[i*N + j]) {
						// 	fill = !fill;
						// 	continue;
						// }

						if (fill) {
							grid[i*N + j] = true;
						}
					}
				}
				break;
				
			case Marching.Vertical:
				for (int j = 0; j < N; j++) {
					bool fill = grid[j];
					
					for (int i = N-1; i >= 0; i--) {
						if (grid[i*N + j]) {
							fill = !fill;
							continue;
						}

						if (fill) {
							grid[i*N + j] = true;
						}
					}
				}
				break;
		}
	}

	public static List<Vec2> MeshFromStencil (bool[] grid, int N, BoundingBox bounds)
	{
		List<Vec2> nodes = new List<Vec2>();
		
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				if (grid[i*N + j]) {
					Vec2 node = ToCartesianSystem(grid, N, i, j, bounds);
					nodes.Add(node);					
				}
			}
		}

		return nodes;
	}

	public static int GetIndex (bool[] grid, int N, int I, int J)
	{
		int idx = 0;
		for (int i = 0; i <= I; i++) {
			for (int j = 0; j <= J; j++) {
				if (grid[i*N + j]) {
					idx += 1;
				}
			}
		}

		return idx;
	}

	public static Gradient GetGradientAtNode (bool[] grid, int N, int i, int j, BoundingBox bounds, double[] values)
	{
		double x_min = bounds.v_min.X;
		double y_min = bounds.v_min.Y;
		double dx = (bounds.v_max.X - bounds.v_min.X) / (double)N;
		double dy = (bounds.v_max.Y - bounds.v_min.Y) / (double)N;		
		double dr = Math.Sqrt(dx*dx + dy*dy);
		double node_value = values[GetIndex(grid, N, i,j)];
		Gradient grad = new Gradient();

		if (grid[i*N + j-1]) {
			grad.CL = (node_value - values[GetIndex(grid, N, i, j-1)]) / dx;
		}
		if (grid[i*N + j+1]) {
			grad.CR = (node_value - values[GetIndex(grid, N, i, j+1)]) / dx;			
		}

		if (grid[(i-1)*N + j-1]) {
			grad.LL = (node_value - values[GetIndex(grid, N, i-1, j-1)]) / dr;
		}
		if (grid[(i-1)*N + j]) {
			grad.LC = (node_value - values[GetIndex(grid, N, i-1, j)]) / dy;
		}
		if (grid[(i-1)*N + j+1]) {
			grad.LR = (node_value - values[GetIndex(grid, N, i-1, j+1)]) / dr;
		}

		if (grid[(i+1)*N + j-1]) {
			grad.UL = (node_value - values[GetIndex(grid, N, i+1, j-1)]) / dr;
		}
		if (grid[(i+1)*N + j]) {
			grad.UC = (node_value - values[GetIndex(grid, N, i+1, j)]) / dy;
		}
		if (grid[(i+1)*N + j+1]) {
			grad.UR = (node_value - values[GetIndex(grid, N, i+1, j+1)]) / dr;
		}

		return grad;
	}

	// ############################################################################
	// Ignore the rest of the function. They are the previous implementation
	// ############################################################################
	public static bool[] GetGrid (Vec2[] vertices, int N)
	{
		var bounds = GetBounds(vertices);
		var grid = new bool[(N+1) * (N+1)];

		var x_min = bounds.v_min.X;
		var y_min = bounds.v_min.Y;
		var dx = (bounds.v_max.X - x_min) / (double)N;
		var dy = (bounds.v_max.Y - y_min) / (double)N;		

		// define the boundaries of on the grid
		for (int i = 0; i < vertices.Length-1; i++) {
			Vec2 a = vertices[i+0];
			Vec2 b = vertices[i+1];
			Vec2 c = Geo.Center(a,b);
			Vec2 n = Geo.Normal(a,b);
			Vec2 l = new Vec2{X = c.X + n.X, Y = c.Y + n.Y};

			int ai = (int)Math.Round((a.Y - y_min) / dy, 0);
			int bi = (int)Math.Round((b.Y - y_min) / dy, 0);
			int ci = (int)Math.Round((c.Y - y_min) / dy, 0);
			int ni = (int)Math.Round((l.Y - y_min) / dy, 0);

			int aj = (int)Math.Round((a.X - x_min) / dx, 0);
			int bj = (int)Math.Round((b.X - x_min) / dx, 0);
			int cj = (int)Math.Round((c.X - x_min) / dx, 0);
			int nj = (int)Math.Round((l.X - x_min) / dx, 0);

			// Console.WriteLine("i: {0}, j:{1}", ci, cj);

			int idx_a = ai*N + aj;
			int idx_b = bi*N + bj;
			int idx_c = ci*N + cj;
			int idx_n = ni*N + nj;

			grid[idx_a] = true;
			grid[idx_b] = true;
			grid[idx_c] = true;
			grid[idx_n] = true;
		}
		
		// fill internal points inside the grid
		for (int i = 0; i < N; i++) {			
			int lhs = 0;
			while (lhs < N && !grid[i*N + lhs]) lhs += 1;

			int rhs = N-1;
			while (rhs > lhs && !grid[i*N + rhs]) rhs -= 1;			

			for (int j = lhs+1; j <= rhs-1; j++) {
				 grid[i*N + j] = true; 					
			}

			// for (int j = 1; j < N/2; j++) {
			// 	if (grid[i*N + j-1]) grid[i*N + j] = true;
			// }
			// for (int j = N-2; j >= N/2; j--) {
			// 	if (grid[i*N + j+1]) grid[i*N + j] = true;
			// }
		}		

		return grid;
	}

	public static Vec2[] ReadFromFile (string path)
	{
		var lines = File.ReadAllLines(path);
		var len   = lines.Length;
		var vertices = new List<Vec2>(len);

		for (int i = 0; i < len; i++) {
			var values = lines[i].Split(',');
			if (values.Length >= 2) {
				var x = Double.Parse(values[1]);
				var y = Double.Parse(values[2]);				
				vertices.Add(new Vec2(x,y));
			}
		}

		return vertices.ToArray();
	}
}

public static class Writer
{
	public static void WriteBoundaryStencil (string path, Vec2[] boundary_nodes, int N)
	{
		var grid = Geo.StencilFromBoundaryNodes(boundary_nodes, N);
		using var fs = System.IO.File.CreateText(path);		

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				if (grid[i*N + j]) {
					fs.Write($"({i},{j})");
				}
				else {
					fs.Write(" ");					
				}		
			}
			fs.Write("\n");
		}
	}

	public static void WriteVertices (string path, bool[] grid, int N, BoundingBox bounds)
	{
		var fs = System.IO.File.CreateText(path);
		var x_min = bounds.v_min.X;
		var y_min = bounds.v_min.Y;
		var dx = (bounds.v_max.X - bounds.v_min.X) / (double)N;
		var dy = (bounds.v_max.Y - bounds.v_min.Y) / (double)N;		
		
		for (int i = 0; i < N; i++) {			
			for (int j = 1; j < N; j++) {
				if (grid[i*N + j] && grid[i*N + j-1]) {
					var x1 = (j-1) * dx + x_min;
					var y1 = (i) * dy + y_min;
					var p1 = new Vec2(x1,y1);
					
					var x2 = (j-0) * dx + x_min;
					var y2 = (i) * dy + y_min;
					var p2 = new Vec2(x2,y2);

					fs.WriteLine($"{p1.X}  {p1.Y}");
					fs.WriteLine($"{p2.X}  {p2.Y}");
					fs.WriteLine("\n");
				}
			}
		}
		
		for (int i = 1; i < N; i++) {
			for (int j = 0; j < N; j++) {			
				if (grid[i*N + j] && grid[(i-1)*N + j]) {
					var x1 = (j) * dx + x_min;
					var y1 = (i-1) * dy + y_min;
					var p1 = new Vec2(x1,y1);
					
					var x2 = (j) * dx + x_min;
					var y2 = (i-0) * dy + y_min;
					var p2 = new Vec2(x2,y2);

					fs.WriteLine($"{p1.X}  {p1.Y}");
					fs.WriteLine($"{p2.X}  {p2.Y}");
					fs.WriteLine("\n");
				}
			}
		}

		fs.Close();
	}
}

class QuadTree
{
	public static bool Traverse (Vec2[] vertices, out int index)
	{
		index = 0;
		
		return false;
	}
}


public class Mesh
{
	private int N;
	private BoundingBox bounds;
	public bool[] grid;
	public List<Vec2> vertices;
	public List<uint> indices;

	public static Mesh FromVertices (Vec2[] vertices, int N)
	{
		var points = new List<Vec2>(N*N);
		var bounds = Geo.GetBounds(vertices);
		var grid   = Geo.GetGrid(vertices, N);

		var x_min = bounds.v_min.X;
		var y_min = bounds.v_min.Y;
		var dx = (bounds.v_max.X - x_min) / (double)N;
		var dy = (bounds.v_max.Y - y_min) / (double)N;		
		var t_filled = 0;
		
		for (int i = 0; i < N; i++) {			
			var y = i * dy + bounds.v_min.Y;

			for (int j = 0; j < N; j++) {
				var x = j * dx + bounds.v_min.X;

				if (grid[i*N + j]) {
					var p = new Vec2(x,y);
					t_filled += 1;
					points.Add(p);
				}
			}
		}
		Console.WriteLine("t_filled: {0}", t_filled);

		return new Mesh{
			N = N,
			bounds = bounds,
			grid = grid,
			vertices = points,
			indices = new(),
		};
	}

	// static bool IsInternalNode (Vec2 p, Vec2[] vertices, BoundingBox bounds)
	// {
	// 	var is_internal_node = false;
	// 	var v_min = bounds.v_min;
	// 	var v_max = bounds.v_max;
	// 	if (p.X < v_min.X || p.X > v_max.X) return false;
	// 	if (p.Y < v_min.Y || p.Y > v_max.Y) return false;

	// 	for (int i = 1; i < vertices.Length; i++) {
	// 		var a = vertices[i-1];
	// 		var b = vertices[i-0];
	// 		var t = Geo.Tangent(a, b);
	// 		var n = Geo.Normal(a, b);

	// 		if (n.X > 0 && n.Y > 0) {
	// 			if (p.X > a.X && p.Y < b.Y) is_internal_node = true;			
	// 		}

	// 		if (n.X < 0 && n.Y > 0) {
	// 			if (p.X < b.X && p.Y > a.Y) is_internal_node = true;
	// 		}

	// 		if (n.X > 0 && n.Y < 0) {
	// 			if (p.X > a.X && p.Y < b.Y) is_internal_node = true;
	// 		}

	// 		if (n.X < 0 && n.Y < 0) {
	// 			if (p.X > b.X && p.Y < a.Y) is_internal_node = true;
	// 		}
	// 	}

	// 	return is_internal_node;
	// }

	/// <summary>N is the number of internal offsets/ levels to create layers of mesh</summary>
	// public static Mesh FromFile(string path, int N)
	// {
	// 	var lines = File.ReadAllLines(path);
	// 	var len   = lines.Length;
	// 	var mesh  = new Mesh();

	// 	for (int i = 0; i < len; i++) {
	// 		var values = lines[i].Split(',');
	// 		if (values.Length >= 2) {
	// 			var x = Double.Parse(values[1]);
	// 			var y = Double.Parse(values[2]);				
	// 			mesh.vertices.Add(new Vec2(x,y));
	// 		}
	// 	}

	// 	len = mesh.vertices.Count;
	// 	double x_min = mesh.vertices[0].X;
	// 	double x_max = mesh.vertices[0].X;
	// 	double y_min = mesh.vertices[0].Y;
	// 	double y_max = mesh.vertices[0].Y;

	// 	for (int i = 0; i < len; i++) {
	// 		x_min = x_min > mesh.vertices[i].X ? mesh.vertices[i].X : x_min;
	// 		x_max = x_max < mesh.vertices[i].X ? mesh.vertices[i].X : x_max;
	// 		y_min = y_min > mesh.vertices[i].Y ? mesh.vertices[i].Y : y_min;
	// 		y_max = y_max < mesh.vertices[i].Y ? mesh.vertices[i].Y : y_max;
	// 	}

	// 	// normalize to [0,1]
	// 	for (int i = 0; i < len; i++) {
	// 		var x = (mesh.vertices[i].X - x_min) / (x_max - x_min);
	// 		var y = (mesh.vertices[i].Y - y_min) / (y_max - y_min);
	// 		mesh.vertices[i] = new Vec2(x,y);
	// 	}
		
	// 	var n_dx = (x_max - x_min) / (2.0 * N);  // shape should be convex, use half of N, half of 'diameter', i.e like radius 
	// 	var n_dy = (y_max - y_min) / (2.0 * N);

	// 	for (int n = 0; n < N; n++) {
	// 		int count = mesh.vertices.Count;
	// 		for (int i = 1; i < len; i++) {
	// 			var p1 = mesh.vertices[count - len + i-1];
	// 			var p2 = mesh.vertices[count - len + i+0];
	// 			// var p1 = mesh.vertices[i + count - len + 0];
	// 			// var p2 = mesh.vertices[i + count - len + 1];
	// 			var pn = Geo.Normal(p1, p2);    // Calculate Normal times dn (some factor)
	// 			// pn.X *= n_dx;
	// 			// pn.Y *= n_dy;
	// 			var pc = Geo.Center(p1, p2);
	// 			pn.X += pc.X;
	// 			pn.Y += pc.Y;
			
	// 			mesh.vertices.Add(pn);

	// 			mesh.indices.Add((uint)(i + count));
	// 			mesh.indices.Add((uint)(i + count - len + n + 0));
	// 			mesh.indices.Add((uint)(i + count - len + n + 1));
	// 		}			
	// 	}

	// 	return mesh;
	// }

	// totally reduntant way to convert to string for C#
	// in C++ do sprintf 
	static string to_string<T> (T obj) where T : struct
	{
		return obj.ToString();
	}

	public void WriteAll(StreamWriter fs)
	{
		for (int i = 0; i < vertices.Count; i++) {
			string line = to_string(vertices[i].X) + "  " + to_string(vertices[i].Y) + "\n";
			fs.Write(line);
		}

		fs.Write("\n\n");

		for (int i = 0; i < indices.Count; i += 3) {
			string line = to_string(indices[i+0]) + "  " + to_string(indices[i+1]) + "  " + to_string(indices[i+2]) + "\n";
			fs.Write(line);
		}
	}
	
	public void WriteVertices(StreamWriter fs)
	{
		// for (int i = 1; i < vertices.Count; i++) {
		// 	string line0 = to_string(vertices[i-1].X) + "  " + to_string(vertices[i-1].Y) + "\n";
		// 	string line1 = to_string(vertices[i-0].X) + "  " + to_string(vertices[i-0].Y) + "\n";
		// 	fs.Write(line0);
		// 	fs.Write(line1);
		// 	fs.Write("\n");
		// }

		var x_min = bounds.v_min.X;
		var y_min = bounds.v_min.Y;
		var dx = (bounds.v_max.X - bounds.v_min.X) / (double)N;
		var dy = (bounds.v_max.Y - bounds.v_min.Y) / (double)N;		
		
		for (int i = 0; i < N; i++) {			
			for (int j = 1; j < N; j++) {
				if (grid[i*N + j] && grid[i*N + j-1]) {
					var x1 = (j-1) * dx + x_min;
					var y1 = (i) * dy + y_min;
					var p1 = new Vec2(x1,y1);
					
					var x2 = (j-0) * dx + x_min;
					var y2 = (i) * dy + y_min;
					var p2 = new Vec2(x2,y2);

					fs.WriteLine($"{p1.X}  {p1.Y}");
					fs.WriteLine($"{p2.X}  {p2.Y}");
					fs.WriteLine("\n");
				}
			}
		}
		for (int i = 1; i < N; i++) {
			for (int j = 0; j < N; j++) {			
				if (grid[i*N + j] && grid[(i-1)*N + j]) {
					var x1 = (j) * dx + x_min;
					var y1 = (i-1) * dy + y_min;
					var p1 = new Vec2(x1,y1);
					
					var x2 = (j) * dx + x_min;
					var y2 = (i-0) * dy + y_min;
					var p2 = new Vec2(x2,y2);

					fs.WriteLine($"{p1.X}  {p1.Y}");
					fs.WriteLine($"{p2.X}  {p2.Y}");
					fs.WriteLine("\n");
				}
			}
		}
	}

	public void WriteIndices(StreamWriter fs)
	{
		for (int i = 0; i < indices.Count; i += 3) {
			string line = to_string(indices[i+0]) + "  " + to_string(indices[i+1]) + "  " + to_string(indices[i+2]) + "\n";
			fs.Write(line);
		}
	}
}
