namespace GridGeneration;
using System;

using Fn = Func<double, double>; // alias

public struct Sides
{
	public Vec2[] L;
	public Vec2[] R;
	public Vec2[] U;
	public Vec2[] D;
}

public static class GeoRandomizer
{
	// public static Sides RandomSides (Fn fn_l, Fn fn_r, Fn fn_u, Fn fn_d, int N)
	public static Sides RandomSides (Vec2[] l, Vec2[] r, Vec2[] d, Vec2[] u)
	{
		// Vec2[] l = new Vec2[N];
		// Vec2[] r = new Vec2[N];
		// Vec2[] u = new Vec2[N];
		// Vec2[] d = new Vec2[N];
		Random R = Random.Shared;
		int N = u.Length;
		double x_min = l.Min(v => v.X);
		double x_max = r.Max(v => v.X);
		double y_min = d.Min(v => v.Y);
		double y_max = d.Max(v => v.Y);

		double dx = (x_max - x_min) / N;
		double dy = (y_max - y_min) / N;

		for (int i = 0; i < N; i++) {
			// double x = x_min + i * dx;
			// double y = y_min + i * dy;
			
			// double noise_x = dx * 0.5 - R.NextDouble();   // +- 0.5
			// l[i].X += noise_x;

			// noise_x = dx * 0.5 - R.NextDouble();   // +- 0.5
			// r[i].X += noise_x;

			// double noise_y = dy * 0.5 - R.NextDouble();   // +- 0.5
			// u[i].Y += noise_y;

			// noise_y = dy * 0.5 - R.NextDouble();   // +- 0.5
			// d[i].Y += noise_y;
		}

		d[0] = l[0]; 
		d[N-1] = r[0];

		u[0] = l[N-1];
		u[N-1] = r[N-1];

		// fix edge nodes
		

		return new Sides{
			L = l,
			R = r,
			U = u,
			D = d,
		};
	}

	public static Vec2[] SidesToBoundaryNodes (Sides sides, int N)
	{
		Vec2[] nodes = new Vec2[N + N + N + N];
		int k = 0;

		// for (int i = 0; i < N; i++) {
		// 	nodes[i + 0*N] = sides.L[i];
		// 	nodes[i + 3*N] = sides.D[i];
		// }
		// for (int j = 0, i = N-1; i >= 0; i--, j++) {
		// 	nodes[j + 1*N] = sides.R[i];
		// 	nodes[j + 2*N] = sides.U[i];			
		// }
		
		// for (int i = 0; i < N; i++, k++) {
		// 	nodes[k] = sides.L[i];
		// }
		// for (int i = 0; i < N; i++, k++) {
		// 	nodes[k] = sides.D[i];
		// }
		// for (int i = N-1; i > 0; i--, k++) {
		// 	nodes[k] = sides.R[i];
		// }
		// for (int i = N-1; i > 0; i--, k++) {
		// 	nodes[k] = sides.U[i];
		// }

		for (int i = 0; i < N; i++, k++) {
			nodes[k] = sides.L[i];
		}
		for (int i = 0; i < N; i++, k++) {
			nodes[k] = sides.D[i];
		}
		for (int i = 0; i < N; i++, k++) {
			nodes[k] = sides.R[i];
		}
		for (int i = 0; i < N; i++, k++) {
			nodes[k] = sides.U[i];
		}

		return nodes;
	}

	public static Vec2[] BoundsFromMesh (Vec2[,] mesh, int Imax, int Jmax)
	{
		// Vec2[] boundarary_nodes = new Vec2[2*Imax + 2*Jmax];
		List<Vec2> nodes = new List<Vec2>(1000);

		for (int i = 0; i < Imax; i++) {
			nodes.Add(mesh[i,0]);
		}
		for (int j = 0; j < Jmax; j++) {
			nodes.Add(mesh[0,j]);	
		}
		for (int i = Imax-1; i >= 0; i--) {
			nodes.Add(mesh[i,Jmax-1-1]);
		}
		for (int j = Jmax-1; j >= 0; j--) {
			nodes.Add(mesh[Imax-1-1,j]);	
		}

		return nodes.ToArray();
	}
	
	public static void WriteBoundaryNodes (Sides sides, string path, int N)
	{
		using var fs = System.IO.File.CreateText(path);

		for (int i = 0; i < N; i++) {
			Vec2 p = sides.L[i];
			fs.WriteLine(p.X.ToString() + "  " + p.Y.ToString());
		}
		
		fs.Write("\n");
		for (int i = 0; i < N; i++) {
			Vec2 p = sides.R[i];
			fs.WriteLine(p.X.ToString() + "  " + p.Y.ToString());
		}
		
		fs.Write("\n");
		for (int i = 0; i < N; i++) {
			Vec2 p = sides.U[i];
			fs.WriteLine(p.X.ToString() + "  " + p.Y.ToString());
		}
		
		fs.Write("\n");
		for (int i = 0; i < N; i++) {
			Vec2 p = sides.D[i];
			fs.WriteLine(p.X.ToString() + "  " + p.Y.ToString());
		}
		fs.Write("\n");
		
		fs.Close();
	}

	public static Vec2[,] MeshFromRandomizer (Sides sides, int Imax, int Jmax)
	{
		Vec2[,] mesh = new Vec2[Imax,Jmax];

		// copy boundaries
		for (int i = 0; i < Imax; i++) {
			mesh[i,0] = sides.L[i];
			mesh[i,Jmax-1] = sides.R[i];
		}
		for (int j = 0; j < Jmax; j++) {
			mesh[0,j] = sides.D[j];
			mesh[Imax-1,j] = sides.U[j];
		}
		
		// for (int i = 2; i < Imax - 1; i++) {
		// 	Vec2 dξ = (sides.R[i] - sides.L[i]) / (double)(Jmax-1);

		// 	for (int j = 2; j < Jmax - 1; j++) {
		// 		Vec2 dη = (sides.U[j] - sides.D[j]) / (double)(Imax-1);

		// 		mesh[1,1] = mesh[0,0] + dξ + dη;
				

		// 		mesh[i,j] = mesh[i-1,j-1] + dξ + dη;
		// 		mesh[i,j-1].Y = mesh[i,j].Y - dη.Y;
		// 		mesh[i,j+1].Y = mesh[i,j].Y + dη.Y;
		// 	}
		// }
		
		for (int i = 1; i < Imax; i++) {
			Vec2 dξ = (sides.R[i] - sides.L[i]) / (double)(Jmax-1);

			for (int j = 1; j <= i; j++) {
				Vec2 dη = (sides.U[j] - sides.D[j]) / (double)(Imax-1);

				mesh[i,j].X = mesh[i,j-1].X + dξ.X;
				mesh[i,j].Y = mesh[i-1,j].Y + dη.Y; 								
			}
		}
		for (int i = Imax-2; i >= 0; i--) {
			Vec2 dξ = (sides.L[i] - sides.R[i]) / (double)(Jmax-1);
			
			for (int j = Jmax-2; j >= i; j--) {
				Vec2 dη = (sides.D[j] - sides.U[j]) / (double)(Imax-1);

				mesh[i,j].X = mesh[i,j+1].X + dξ.X;
				mesh[i,j].Y = mesh[i+1,j].Y + dη.Y; 				
			}
		}

		return mesh;
	}

	public static Gradient GetGradient (Vec2[,] mesh, double[,] signal, int i, int j)
	{
		Vec2 p = mesh[i,j];
		double s = signal[i,j];

		return new Gradient{
			LL = (s - signal[i-1,j-1]) / (p - mesh[i-1,j-1]).Length(),
			LC = (s - signal[i-1,j+0]) / (p - mesh[i-1,j+0]).Length(),
			LR = (s - signal[i-1,j+1]) / (p - mesh[i-1,j+1]).Length(),
			
			UL = (s - signal[i+1,j-1]) / (p - mesh[i-1,j-1]).Length(),
			UC = (s - signal[i+1,j+0]) / (p - mesh[i-1,j+0]).Length(),
			UR = (s - signal[i+1,j+1]) / (p - mesh[i-1,j+1]).Length(),

			CL = (s - signal[i,j-1]) / (p - mesh[i,j-1]).Length(),
			CR = (s - signal[i,j+1]) / (p - mesh[i,j+1]).Length(),
		};
	}

	public static void WriteCorrectGrid (string path, Vec2[,] mesh, int Imax, int Jmax)
	{
		var fs = System.IO.File.CreateText(path);

		for (int i = 0; i < Imax; i++) {			
			for (int j = 1; j < Jmax; j++) {
				Vec2 p1 = mesh[i,j-1];					
				Vec2 p2 = mesh[i,j-0];

				fs.WriteLine($"{p1.X}  {p1.Y}");
				fs.WriteLine($"{p2.X}  {p2.Y}");
				fs.WriteLine("\n");
			}
		}
		
		for (int i = 1; i < Imax; i++) {
			for (int j = 0; j < Jmax; j++) {			
				Vec2 p1 = mesh[i-1,j];					
				Vec2 p2 = mesh[i-0,j];

				fs.WriteLine($"{p1.X}  {p1.Y}");
				fs.WriteLine($"{p2.X}  {p2.Y}");
				fs.WriteLine("\n");
			}
		}

		fs.Close();
	}

	public static void CheckUnitialized (Vec2[,] mesh, int Imax, int Jmax)
	{
		for (int i = 0; i < Imax; i++) {
			for (int j = 0; j < Jmax; j++) {
				if (mesh[i,j].X == 0 && mesh[i,j].Y == 0) {
					Console.WriteLine("unitialized: ({0}, {1})", i, j);
				}
			}
		}
	}
}
