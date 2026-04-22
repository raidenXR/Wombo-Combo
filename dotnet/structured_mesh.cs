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

public struct Grid
{
	public bool[] stencil_buffer;
	public int[] indices;
	public int N;
	public BoundingBox bounds;

	public Grid(Vec2[] boundary_nodes, int N) {	
		stencil_buffer = new bool[N*N];
		indices = new int[N];
		this.N = N;
		bounds = Geo.GetBounds(boundary_nodes);	
	}
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
	public static Vec2 ToCartesianSystem (Grid grid, int i, int j)
	{
		double x_min = grid.bounds.v_min.X;
		double y_min = grid.bounds.v_min.Y;
		double dx = (grid.bounds.v_max.X - grid.bounds.v_min.X) / (double)grid.N;
		double dy = (grid.bounds.v_max.Y - grid.bounds.v_min.Y) / (double)grid.N;		

		return new Vec2{
			X = (double)j * dx + x_min,
			Y = (double)i * dy + y_min,
		};		
	}

	/// <summary> transforms to i,j (Stencilsystem) from x,y (Cartesian system) </summary> 
	public static Index ToStencilSystem (int N, Vec2 p, BoundingBox bounds)
	{
		double x_min = bounds.v_min.X;
		double x_max = bounds.v_max.X;
		double y_min = bounds.v_min.Y;
		double y_max = bounds.v_max.Y;
		// double dx = (x_max - x_min) / (double)N;
		// double dy = (y_max - y_min) / (double)N;		
		
		int i = (int)Math.Round((N-1)*(p.Y - y_min) / (y_max - y_min), 0); // normalize to [0..N-1]
		int j = (int)Math.Round((N-1)*(p.X - x_min) / (x_max - x_min), 0); // normalize to [0..N-1]

		// return new Index{
		// 	I = Math.Clamp(i, 0, N-1),
		// 	J = Math.Clamp(j, 0, N-1),
		// };		
		return new Index{
			I = i,
			J = j,
		};		
	}

	/// <summary> -- RECURSIVE FUNCTION --
	/// stencil method-algorithm for setting to 1 (true) the elements that are part of the boundary
	/// of the domain </summary>
	public static void AssignStencilElement (Grid grid, Vec2 a, Vec2 b)
	{
		var index_a = ToStencilSystem(grid.N, a, grid.bounds);
		var index_b = ToStencilSystem(grid.N, b, grid.bounds);

		// subdived (recursively) if nessacary
		if (Math.Abs(index_a.I - index_b.I) > 1 || Math.Abs(index_a.J - index_b.J) > 1) {
			AssignStencilElement(grid, a, Center(a,b));
			AssignStencilElement(grid, Center(a,b), b);
		}		

		grid.stencil_buffer[index_a.I * grid.N + index_a.J] = true;
		grid.stencil_buffer[index_b.I * grid.N + index_b.J] = true;
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

	public static Grid StencilFromBoundaryNodes (Vec2[] boundary_nodes, int N)
	{
		int len = boundary_nodes.Length;
		Grid grid = new Grid(boundary_nodes, N);
		
		for (int i = 0; i < len-1; i++) {
			Vec2 a = boundary_nodes[i+0];
			Vec2 b = boundary_nodes[i+1];
			AssignStencilElement(grid, a, b);
			
			Vec2 c = Center(a,b);
			Index ci = ToStencilSystem(N, c, grid.bounds);
			grid.stencil_buffer[ci.I*N + ci.J] = true;
		}		

		return grid;
	}

	static int MeasureMarchingRow (bool[] stencil_buffer, int N, int I)
	{
		int n = 0;
		for (int j = 0; j < N; j++) {
			if (stencil_buffer[I*N + j]) {
				while (stencil_buffer[I*N + j]) j++;  // advance
				n += 1;
			}
		}

		return n;
	}

	public static void FillStencil (Grid grid)
	{
		for (int i = 0; i < grid.N; i++) {
			int lhs = 0;
			int rhs = grid.N-1;

			while (!grid.stencil_buffer[i*grid.N + lhs]) lhs++;  // advance 
			while (!grid.stencil_buffer[i*grid.N + rhs]) rhs--;  // advance

			while (grid.stencil_buffer[i*grid.N + lhs]) lhs++;   // advance
			while (grid.stencil_buffer[i*grid.N + rhs]) rhs--;   // advance

			int a = Math.Min(lhs, rhs);
			int b = Math.Max(lhs, rhs);
			bool fill = true;
			
			switch (MeasureMarchingRow(grid.stencil_buffer, grid.N, i)) {
				case 2:
				case 3:
					for (int j = a; j <= b; j++) {
						grid.stencil_buffer[i*grid.N + j] = true;
					}
					break;

				case 4: 
					for (int j = a; j <= b; j++) {
						if (grid.stencil_buffer[i*grid.N + j]) {
							while (grid.stencil_buffer[i*grid.N + j]) j++;  // advance
							j--;
							fill = !fill;
							continue;
						}

						if (fill) {
							grid.stencil_buffer[i*grid.N + j] = true;
						}
					}
				break;
				default: break;
			}
		}

		int offset = 0;
		for (int i = 0; i < grid.N; i++) {
			grid.indices[i] = offset;
			
			for (int j = 0; j < grid.N; j++) {
				if (grid.stencil_buffer[i*grid.N + j]) {
					offset += 1;
				}
			}
		}
	}

	static void FillStencil_element_Y (Grid grid, Vec2 a, Vec2 b)
	{		
		var index_a = ToStencilSystem(grid.N, a, grid.bounds);
		var index_b = ToStencilSystem(grid.N, b, grid.bounds);

		// subdived (recursively) if nessacary
		if (Math.Abs(index_a.I - index_b.I) > 1 || Math.Abs(index_a.J - index_b.J) > 1) {
			FillStencil_element_Y(grid, a, Center(a,b));
			FillStencil_element_Y(grid, Center(a,b), b);
		}		

		int N = grid.N;		
		Vec2 c = Center(a, b);
		Vec2 n = Normal(a, b);

		Index idx = ToStencilSystem(N, a, grid.bounds);
		// Index idx = ToStencilSystem(N, new Vec2(c.X + n.X, c.Y + n.Y), grid.bounds);
		int I = idx.I;
		int J = idx.J;

		if (n.Y >= 0) {
			for (int i = I+1; i < N && !grid.stencil_buffer[i*N + J]; i++) {
				grid.stencil_buffer[i*N + J] = true;					
			}					
		}
		if (n.Y < 0) {
			for (int i = I-1; i >= 0 && !grid.stencil_buffer[i*N + J]; i--) {
				grid.stencil_buffer[i*N + J] = true;					
			}					
		}
		// if (n.X >= 0) {
		// 	for (int j = J+1; j < N && !grid.stencil_buffer[I*N + j]; j++) {
		// 		grid.stencil_buffer[I*N + j] = true;					
		// 	}					
		// }
		// if (n.X < 0) {
		// 	for (int j = J-1; j >= 0 && !grid.stencil_buffer[I*N + j]; j--) {
		// 		grid.stencil_buffer[I*N + j] = true;					
		// 	}					
		// }
	}

	static void FillStencil_element_X (Grid grid, Vec2 a, Vec2 b)
	{		
		var index_a = ToStencilSystem(grid.N, a, grid.bounds);
		var index_b = ToStencilSystem(grid.N, b, grid.bounds);

		// subdived (recursively) if nessacary
		if (Math.Abs(index_a.I - index_b.I) > 1 || Math.Abs(index_a.J - index_b.J) > 1) {
			FillStencil_element_X(grid, a, Center(a,b));
			FillStencil_element_Y(grid, Center(a,b), b);
		}		

		int N = grid.N;		
		Vec2 c = Center(a, b);
		Vec2 n = Normal(a, b);

		Index idx = ToStencilSystem(N, a, grid.bounds);
		// Index idx = ToStencilSystem(N, new Vec2(c.X + n.X, c.Y + n.Y), grid.bounds);
		int I = idx.I;
		int J = idx.J;

		if (n.X >= 0) {
			for (int j = J+1; j < N && !grid.stencil_buffer[I*N + j]; j++) {
				grid.stencil_buffer[I*N + j] = true;					
			}					
		}
		if (n.X < 0) {
			for (int j = J-1; j >= 0 && !grid.stencil_buffer[I*N + j]; j--) {
				grid.stencil_buffer[I*N + j] = true;					
			}					
		}
	}

	public static void FillStencil_2 (Grid grid, Vec2[] boundary_nodes)
	{
		int N = grid.N;		
		
		for (int k = 0; k < N-1; k++) {
			Vec2 a = boundary_nodes[k+0];
			Vec2 b = boundary_nodes[k+1];
			Vec2 c = Center(a, b);
			Vec2 n = Normal(a, b);

			FillStencil_element_Y(grid, a, b);
		}

		for (int k = 0; k < N-1; k++) {
			Vec2 a = boundary_nodes[k+0];
			Vec2 b = boundary_nodes[k+1];
			Vec2 c = Center(a, b);
			Vec2 n = Normal(a, b);
			FillStencil_element_X(grid, a, b);

		}

		int offset = 0;
		for (int i = 0; i < grid.N; i++) {
			grid.indices[i] = offset;
			
			for (int j = 0; j < grid.N; j++) {
				if (grid.stencil_buffer[i*grid.N + j]) {
					offset += 1;
				}
			}
		}
	}

	public static void Reverse (Grid grid)
	{
		int N = grid.N;
		
		for (int i = 0; i < N; i++)	{
			for (int j = 0; j < N; j++) {
				grid.stencil_buffer[i*N + j] = !grid.stencil_buffer[i*N + j];
			}
		}
	}
	
	public static List<Vec2> MeshFromStencil (Grid grid)
	{
		List<Vec2> nodes = new List<Vec2>();
		
		for (int i = 0; i < grid.N; i++) {
			for (int j = 0; j < grid.N; j++) {
				if (grid.stencil_buffer[i*grid.N + j]) {
					Vec2 node = ToCartesianSystem(grid, i, j);
					nodes.Add(node);					
				}
			}
		}

		return nodes;
	}

	public static int GetIndex (Grid grid, int I, int J)
	{
		// int idx = grid.indices[I];
		// for (int j = 0; j <= J; j++) {
		// 	idx += (grid.stencil_buffer[I*grid.N + j]) ? 1 : 0;
		// }

		// return idx;
		
		int idx = 0;
		for (int i = 0; i <= I; i++) {
			for (int j = 0; j <= J; j++) {
				idx += (grid.stencil_buffer[I*grid.N + j]) ? 1 : 0;
			}			
		}

		return idx;
	}

	public static Gradient GetGradientAtNode (Grid grid, int i, int j, double[] values)
	{
		// if (i < 0 || i > grid.N-1) throw new ArgumentException($"{i} is not a valid Index input");
		// if (j < 0 || j > grid.N-1) throw new ArgumentException($"{j} is not a valid Index input");

		int N = grid.N;
		if (i < 0 || i > N-1 || j < 0 || j > N-1) {
			return new Gradient();
		}
		
		BoundingBox bounds = grid.bounds;
		double x_min = bounds.v_min.X;
		double y_min = bounds.v_min.Y;
		double dx = (bounds.v_max.X - bounds.v_min.X) / (double)N;
		double dy = (bounds.v_max.Y - bounds.v_min.Y) / (double)N;		
		double dr = Math.Sqrt(dx*dx + dy*dy);
		double node_value = values[GetIndex(grid, i,j)];
		Gradient grad = new Gradient();

		if (j > 0 && grid.stencil_buffer[i*N + j-1]) {
			grad.CL = (node_value - values[GetIndex(grid, i, j-1)]) / dx;
		}
		if (j < N-1 && grid.stencil_buffer[i*N + j+1]) {
			grad.CR = (node_value - values[GetIndex(grid, i, j+1)]) / dx;			
		}

		if (j > 0 && i > 0 && grid.stencil_buffer[(i-1)*N + j-1]) {
			grad.LL = (node_value - values[GetIndex(grid, i-1, j-1)]) / dr;
		}
		if (i > 0 && grid.stencil_buffer[(i-1)*N + j]) {
			grad.LC = (node_value - values[GetIndex(grid, i-1, j)]) / dy;
		}
		if (j < N-1 && i > 0 && grid.stencil_buffer[(i-1)*N + j+1]) {
			grad.LR = (node_value - values[GetIndex(grid, i-1, j+1)]) / dr;
		}

		if (j > 0 && i < N-1 && grid.stencil_buffer[(i+1)*N + j-1]) {
			grad.UL = (node_value - values[GetIndex(grid, i+1, j-1)]) / dr;
		}
		if (i < N-1 && grid.stencil_buffer[(i+1)*N + j]) {
			grad.UC = (node_value - values[GetIndex(grid, i+1, j)]) / dy;
		}
		if (j < N-1 && i < N-1 && grid.stencil_buffer[(i+1)*N + j+1]) {
			grad.UR = (node_value - values[GetIndex(grid, i+1, j+1)]) / dr;
		}

		return grad;
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

	public static bool[] ParseStencilBuffer (string path, char delimeter, int N)
	{
		var lines = File.ReadAllLines(path);
		var stencil_buffer = new bool[N*N];

		for (int i = 0; i < N; i++) {
			if (lines[N-i-1].Length > 1) {  // the text file is upside-down !!
				for (int j = 0; j < N; j++) {
					stencil_buffer[i*N + j] = lines[N-i-1][j] == delimeter ? true : false;
				}
			}	
		}

		return stencil_buffer;
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

	public static Vec2 RandomAntennaPosition (List<Vec2> nodes)
	{
		int a = Random.Shared.Next(nodes.Count);
		int b = a;
		
		while (a == b) {
			b = Random.Shared.Next(nodes.Count);
		}
		
		return Geo.Center(nodes[a], nodes[b]);	
	}

	public static Vec2[] InitialAntennasPositions (List<Vec2> nodes, int n)
	{
		Vec2[] positions = new Vec2[n];
		
		for (int i = 0; i < n; i++) {			
			positions[i] = RandomAntennaPosition(nodes);	
		}

		return positions;
	}
	
	public static int GetClosestElement (Vec2 p, Vec2[] boundary_nodes)
	{
		int idx = -1;
		int K = boundary_nodes.Length;
		double r_min = Double.MaxValue;

		for (int k = 0; k < K-1; k++) {
			Vec2 a = boundary_nodes[k+0];
			Vec2 b = boundary_nodes[k+1];
			Vec2 c = Geo.Center(a, b);
			
			double r = Math.Sqrt((c.X - p.X)*(c.X - p.X) + (c.Y - p.Y)*(c.Y - p.Y));
			
			if (r_min > r) {
				r_min = r;
				idx = k;
			}			
		}

		return idx;
	}

	public static int GetClosestPoint (Vec2 p, Vec2[] positions)
	{
		int idx = 0;
		int K = positions.Length;
		double r_min = Math.Sqrt((p.X - positions[0].X) * (p.X - positions[0].X) + (p.Y - positions[0].Y) * (p.Y - positions[0].Y));

		for (int k = 0; k < K; k++) {			
			double r = Math.Sqrt((p.X - positions[k].X) * (p.X - positions[k].X) + (p.Y - positions[k].Y) * (p.Y - positions[k].Y));
			
			if (r_min > r) {
				r_min = r;
				idx = k;
			}			
		}

		return idx;		
	}


	/// <summary> ensure that antennas are INSIDE the domain !!! </summary>
	public static void ConstrainAntennasPositions (Grid grid, Vec2[] antennas_pos, Vec2[] boundary_nodes)
	{
		int N = grid.N;
		BoundingBox bounds = grid.bounds;
		double dx = (bounds.v_max.X - bounds.v_min.X) / (double)N;
		double dy = (bounds.v_max.Y - bounds.v_min.Y) / (double)N;		

		for (int k = 0; k < antennas_pos.Length; k++) {
			Index idx = Geo.ToStencilSystem(N, antennas_pos[k], bounds);
			int i = idx.I;
			int j = idx.J;

			// ensure that antenna is INSIDE the domain !!!
			if ((i < 0 || i > N-1) || (j < 0 || j > N-1)) {
				int i_a = GetClosestElement(antennas_pos[k], boundary_nodes);
				int i_b = i_a + 1;
				Vec2 t = Geo.Tangent(boundary_nodes[i_a], boundary_nodes[i_b]);
				antennas_pos[k].X = t.X;
				antennas_pos[k].Y = t.Y;
				
				idx = Geo.ToStencilSystem(N, antennas_pos[k], bounds);
				i = idx.I;
				j = idx.J;
			}
		}
		
	}

	static void MinMax (double[] values, ref double min, ref double max)
	{
		min = values[0];
		max = values[0];

		for (int i = 0; i < values.Length; i++) {
			min = values[i] < min ? values[i] : min;
			max = values[i] > max ? values[i] : max;
		}
	}

	/// <summary> </summary>
	/// <param name="signal_intensity"> the array that stores the intensities on each node of the mesh </param>
	/// <param name="I_opt"> the intesity value used for the optimization, to compare against </param>
	public static void MoveAntennas (Grid grid, Vec2[] antennas_pos, Vec2[] boundary_nodes, double[] signal_intensity)
	{
		int N = grid.N;
		BoundingBox bounds = grid.bounds;
		double dx = (bounds.v_max.X - bounds.v_min.X) / (double)N;
		double dy = (bounds.v_max.Y - bounds.v_min.Y) / (double)N;		

		double min = 0, max = 0;
		MinMax(signal_intensity, ref min, ref max);

		for (int k = 0; k < antennas_pos.Length; k++) {
			GradDirection dir = GradDirection.Self;
			Index idx = Geo.ToStencilSystem(N, antennas_pos[k], bounds);
			int i = idx.I;
			int j = idx.J;

			Gradient grad = Geo.GetGradientAtNode(grid, i, j, signal_intensity);				

			if (signal_intensity[Geo.GetIndex(grid, i, j)] > 0.7 * max) {
				dir = grad.SteepestDescend();
			}
			else if (signal_intensity[Geo.GetIndex(grid, i, j)] < 0.2  * max) {
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
		}
		ConstrainAntennasPositions(grid, antennas_pos, boundary_nodes);
	}

	/// <summary> </summary>
	/// <param name="signal_intensity"> the array that stores the intensities on each node of the mesh </param>
	/// <param name="I_opt"> the intesity value used for the optimization, to compare against </param>
	public static void MoveAntennas_2 (Grid grid, Vec2[] antennas_pos, Vec2[] boundary_nodes, double[] signal_intensity)
	{
		int N = grid.N;
		BoundingBox bounds = grid.bounds;
		double dx = (bounds.v_max.X - bounds.v_min.X) / (double)N;
		double dy = (bounds.v_max.Y - bounds.v_min.Y) / (double)N;		

		double min = 0, max = 0;
		MinMax(signal_intensity, ref min, ref max);
		double mean = min + (max - min) /2;

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				if (grid.stencil_buffer[i*N + j]) {					
					Gradient grad = Geo.GetGradientAtNode(grid, i, j, signal_intensity);				
					GradDirection dir = GradDirection.Self;
					Vec2 dv = new Vec2();

					if (signal_intensity[Geo.GetIndex(grid, i, j)] > 0.7 * mean) {
						dir = grad.SteepestAscend();
					}
					else if (signal_intensity[Geo.GetIndex(grid, i, j)] < 0.2  * mean) {
						dir = grad.SteepestDescend();
					}

					switch (dir) {
				        case GradDirection.UL:
				            dv.X -= dx; 
				            dv.Y += dy;
							break;
				        case GradDirection.UC:
				            dv.Y += dy;
							break;
				        case GradDirection.UR: 
				            dv.X += dx;
				            dv.Y += dy;
							break;
				        case GradDirection.CL: 
				            dv.X -= dx;
							break;
				        case GradDirection.CR: 
				            dv.X += dx;
							break;
				        case GradDirection.LL: 
				            dv.X -= dx;
				            dv.Y -= dy;
							break;
				        case GradDirection.LC: 
				            dv.Y -= dy;
							break;
				        case GradDirection.LR: 
				            dv.X += dx;
				            dv.Y -= dy;
							break;				
						default:
							break;				
					}

					Vec2 p = Geo.ToCartesianSystem(grid, i, j);
					int k = GetClosestPoint(p, antennas_pos);
					antennas_pos[k].X -= dv.X;
					antennas_pos[k].Y -= dv.Y;
				}
			}
		}

		ConstrainAntennasPositions(grid, antennas_pos, boundary_nodes);		
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
				if (grid.stencil_buffer[i*N + j]) {
					fs.Write($"({i},{j})");
				}
				else {
					fs.Write(" ");					
				}		
			}
			fs.Write("\n");
		}
	}


	public static void WriteBoundaryStencil (string path, Vec2[] boundary_nodes, int N, char delimeter_true, char delimeter_false)
	{
		var grid = Geo.StencilFromBoundaryNodes(boundary_nodes, N);
		using var fs = System.IO.File.CreateText(path);		

		for (int i = N-1; i >= 0; i--) {
			for (int j = 0; j < N; j++) {
				if (grid.stencil_buffer[i*N + j]) {
					fs.Write(delimeter_true);
				}
				else {
					fs.Write(delimeter_false);					
				}		
			}
			fs.Write("\n");
		}
	}
	
	public static void WriteGridStencil (string path, Grid grid, char delimeter_true, char delimeter_false)
	{
		using var fs = System.IO.File.CreateText(path);		
		int N = grid.N;
		
		for (int i = N-1; i >= 0; i--) {
			for (int j = 0; j < N; j++) {
				if (grid.stencil_buffer[i*N + j]) {
					fs.Write(delimeter_true);
				}
				else {
					fs.Write(delimeter_false);					
				}		
			}
			fs.Write("\n");
		}
	}
	
	public static void WriteVertices (string path, Grid grid)
	{
		int N = grid.N;
		BoundingBox bounds = grid.bounds;
		var fs = System.IO.File.CreateText(path);
		double x_min = bounds.v_min.X;
		double y_min = bounds.v_min.Y;
		double dx = (bounds.v_max.X - bounds.v_min.X) / (double)N;
		double dy = (bounds.v_max.Y - bounds.v_min.Y) / (double)N;		
		
		for (int i = 0; i < N; i++) {			
			for (int j = 1; j < N; j++) {
				if (grid.stencil_buffer[i*N + j] && grid.stencil_buffer[i*N + j-1]) {
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
				if (grid.stencil_buffer[i*N + j] && grid.stencil_buffer[(i-1)*N + j]) {
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

