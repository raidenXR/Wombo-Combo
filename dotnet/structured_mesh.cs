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

// public struct Grid
// {
// 	public bool[] grid;
// 	public int N;
// 	public BoundingBox bounds;

// 	public Grid(Vec2[] boundary_nodes, int N) {	
// 		grid = new bool[N*N];
// 		this.N = N;
// 		bounds = Geo.GetBounds(boundary_nodes);	
// 	}
// }

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

		if (j > 0 && grid[i*N + j-1]) {
			grad.CL = (node_value - values[GetIndex(grid, N, i, j-1)]) / dx;
		}
		if (j < N-1 && grid[i*N + j+1]) {
			grad.CR = (node_value - values[GetIndex(grid, N, i, j+1)]) / dx;			
		}

		if (j > 0 && i > 0 && grid[(i-1)*N + j-1]) {
			grad.LL = (node_value - values[GetIndex(grid, N, i-1, j-1)]) / dr;
		}
		if (i > 0 && grid[(i-1)*N + j]) {
			grad.LC = (node_value - values[GetIndex(grid, N, i-1, j)]) / dy;
		}
		if (j < N-1 && i > 0 && grid[(i-1)*N + j+1]) {
			grad.LR = (node_value - values[GetIndex(grid, N, i-1, j+1)]) / dr;
		}

		if (j > 0 && i < N-1 && grid[(i+1)*N + j-1]) {
			grad.UL = (node_value - values[GetIndex(grid, N, i+1, j-1)]) / dr;
		}
		if (i < N-1 && grid[(i+1)*N + j]) {
			grad.UC = (node_value - values[GetIndex(grid, N, i+1, j)]) / dy;
		}
		if (j < N-1 && i < N-1 && grid[(i+1)*N + j+1]) {
			grad.UR = (node_value - values[GetIndex(grid, N, i+1, j+1)]) / dr;
		}

		return grad;
	}

	public static int GetClosestNode (bool[] grid, int N, int i, int j, BoundingBox bounds, Vec2[] nodes)
	{
		int idx = 0;
		int K = nodes.Length;
		double r_min = Double.MaxValue;
		Vec2 p = ToCartesianSystem(grid, N, i, j, bounds);

		for (int k = 0; k < K; k++) {
			Vec2 v = nodes[k];
			double r = Math.Sqrt(v.X*v.X + v.Y*v.Y);
			if (r_min > r) {
				r_min = r;
				idx = k;
			}			
		}

		return idx;
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

public static class Signal
{
	/// <summary> Gauss function for the curve of the signal Intensity  </summary>
	/// <param name="A"> amplitude </param>
	/// <param name="x"> x </param>
	/// <param name="x0"> center -x </param>
	/// <param name="y"> y </param>
	/// <param name="y0"> center -y </param>
	/// <param name="s_x"> variance -x </param>
	/// <param name="s_y"> variance -y </param>
	public static double Gaussian (double A, double x, double x0, double y, double y0, double s_x, double s_y)
	{
		double t1 = (Math.Pow(x-x0, 2)) / (2 * Math.Pow(s_x, 2));
		double t2 = (Math.Pow(y-y0, 2)) / (2 * Math.Pow(s_y, 2));

		return A * Math.Exp(-(t1 + t2));
	}

	public static Vec2[] InitialAntennasPositions (List<Vec2> nodes, int n)
	{
		Vec2[] positions = new Vec2[n];
		
		for (int i = 0; i < n; i++) {
			int a = Random.Shared.Next(nodes.Count);
			int b = a;
			
			while (a == b) {
				b = Random.Shared.Next(nodes.Count);
			}
			
			positions[i] = Geo.Center(nodes[a], nodes[b]);	
		}

		return positions;
	}

	/// <summary> </summary>
	/// <param name="signal_intensity"> the array that stores the intensities on each node of the mesh </param>
	/// <param name="I_opt"> the intesity value used for the optimization, to compare against </param>
	public static void MoveAntennas (bool[] grid, int N, BoundingBox bounds, Vec2[] antennas_pos, double[] signal_intensity)
	{
		double dx = (bounds.v_max.X - bounds.v_min.X) / (double)N;
		double dy = (bounds.v_max.Y - bounds.v_min.Y) / (double)N;		

		for (int k = 0; k < antennas_pos.Length; k++) {
				Index idx = Geo.ToStencilSystem(N, antennas_pos[k], bounds);
				int i = idx.I;
				int j = idx.J;
		// for (int i = 1; i < N-1; i++) {
			// for (int j = 1; j < N-1; j++) {
				// int k =  Geo.GetClosestNode(grid, N, i, j, bounds, antennas_pos);
				Gradient grad = Geo.GetGradientAtNode(grid, N, i, j, bounds, signal_intensity);
				GradDirection dir = GradDirection.Self;

				if (signal_intensity[Geo.GetIndex(grid, N, i, j)] > 8) {
					dir = grad.SteepestDescend();
				}
				else if (signal_intensity[Geo.GetIndex(grid, N, i, j)] < 4) {
					dir = grad.SteepestAscend();
				}
				else {
					dir = GradDirection.Self;
				}

				switch (dir) {
			        case GradDirection.UL:
			            antennas_pos[k].X -= dx; 
			            antennas_pos[k].Y += dy;
						break;
			        case GradDirection.UC:
			            antennas_pos[k].Y += dy;
						break;
			        case GradDirection.UR: 
			            antennas_pos[k].X += dx;
			            antennas_pos[k].Y += dy;
						break;
			        case GradDirection.CL: 
			            antennas_pos[k].X -= dx;
						break;
			        case GradDirection.CR: 
			            antennas_pos[k].X += dx;
						break;
			        case GradDirection.LL: 
			            antennas_pos[k].X -= dx;
			            antennas_pos[k].Y -= dy;
						break;
			        case GradDirection.LC: 
			            antennas_pos[k].Y -= dy;
						break;
			        case GradDirection.LR: 
			            antennas_pos[k].X += dx;
			            antennas_pos[k].Y -= dy;
						break;				
					default:
						break;
				}
			// }
		}
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
		double x_min = bounds.v_min.X;
		double y_min = bounds.v_min.Y;
		double dx = (bounds.v_max.X - bounds.v_min.X) / (double)N;
		double dy = (bounds.v_max.Y - bounds.v_min.Y) / (double)N;		
		
		for (int i = 0; i < N; i++) {			
			for (int j = 1; j < N; j++) {
				if (grid[i*N + j] && grid[i*N + j-1]) {
					double x1 = (j-1) * dx + x_min;
					double y1 = (i) * dy + y_min;
					Vec2 p1 = new Vec2(x1,y1);
					
					double x2 = (j-0) * dx + x_min;
					double y2 = (i) * dy + y_min;
					Vec2 p2 = new Vec2(x2,y2);

					fs.WriteLine($"{p1.X}  {p1.Y}");
					fs.WriteLine($"{p2.X}  {p2.Y}");
					fs.WriteLine("\n");
				}
			}
		}
		
		for (int i = 1; i < N; i++) {
			for (int j = 0; j < N; j++) {			
				if (grid[i*N + j] && grid[(i-1)*N + j]) {
					double x1 = (j) * dx + x_min;
					double y1 = (i-1) * dy + y_min;
					Vec2 p1 = new Vec2(x1,y1);
					
					double x2 = (j) * dx + x_min;
					double y2 = (i-0) * dy + y_min;
					Vec2 p2 = new Vec2(x2,y2);

					fs.WriteLine($"{p1.X}  {p1.Y}");
					fs.WriteLine($"{p2.X}  {p2.Y}");
					fs.WriteLine("\n");
				}
			}
		}

		fs.Close();
	}
}

