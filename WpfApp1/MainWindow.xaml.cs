using HelixToolkit.Wpf;
using MathLib;
using System;
using System.Diagnostics;
using System.Linq;
using System.Windows;
using System.Windows.Media;
using System.Windows.Media.Media3D;

namespace WpfApp1
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        private const double LineDiameter = 0.05;
        private const double LineInnerDiameter = 0;
        private readonly Brush LineFill = Brushes.LightGray;
        private const double VertRadius = 0.1;

        private const double E = 1.0;
        private const double v = 1.0;
        private const double l = E / ((1 + v) * (1 - 2 * v));
        private const double m = E / (2 * (1 + v));

        private readonly double LargeCoefficient = Math.Pow(10, 20);

        public MainWindow()
        {
            InitializeComponent();
         
            const int nx = 2;
            const int ny = 1;
            const int nz = 2;

            const int ax = 3;
            const int ay = 3;
            const int az = 3;

            // initial points array
            var AKT = GenerateAKT(nx, ny, nz, ax, ay, az);
            Render(AKT);

            // TODO: calculate NT
            var NT = GeneratorNP.Generate(globalToAKT, nx, ny, nz);

            // fixed points global coords (ZU) and force info (ZP)
            var ZU = CalculateZU(nx, ny, nz);
            var ZP = CalculateZP(nx, ny, nz);

            var DFIABG = GenerateDFIABG();
            var DPSITE = GenerateDPSITE.Generate();
            
            var (MG, F) = ProcessElements(nx, ny, nz, AKT, NT, DFIABG, DPSITE, ZP);
            FixMG(MG, ZU);

            var U = GaussianElimination(MG, F);
            var result = AddTranslation(AKT, U);
            RenderResult(result);
        }
        
        private void Render(Point3D[] points)
        {
            // TODO: add edges
            GenerateAndAddShapeToCollection(points, new int[] { }, Brushes.Blue, MainViewPort.Children);
        }

        private void RenderResult(Point3D[] points)
        {
            // TODO: add edges
            GenerateAndAddShapeToCollection(points, new int[] { }, Brushes.Red, MainViewPort.Children);
        }

        private static int[,,] globalToAKT;

        private Point3D[] GenerateAKT(int nx, int ny, int nz, double ax, double ay, double az)
        {
            double sx = ax / nx;
            double sy = ay / ny;
            double sz = az / nz;
            int cx = nx + 1;
            int cy = ny + 1;
            int cz = nz + 1;
            int cxyz = cx * cy * cz;                          // 1x1 boxes vertex count
            int ce = (cx * nz + nx * cz) * cy + ny * cx * cz; // edges count

            var points = new Point3D[cxyz + ce];
            globalToAKT = new int[cx * 2 - 1, cy * 2 - 1, cz * 2 - 1];

            int i = 0;
            for (int iz = 0; iz < cz * 2 - 1; ++iz)
            {
                for (int iy = 0; iy < cy * 2 - 1; ++iy)
                {
                    for (int ix = 0; ix < cx * 2 - 1; ++ix)
                    {
                        int count = ix % 2 + iy % 2 + iz % 2;
                        if (count <= 1)
                        {
                            globalToAKT[ix, iy, iz] = i;

                            points[i] = new Point3D(ix * sx / 2.0, iy * sy / 2.0, iz * sz / 2.0);
                            i++;
                        }
                    }
                }
            }

#if DEBUG
            Debug.WriteLine("Initial points (AKT): ");
            foreach (var p in points)
            {
                Debug.WriteLine(p);
            }
#endif

            return points;
        }

        private int[] CalculateZU(int nx, int ny, int nz)
        {
            // basically fixed points are just points on bottom (X-Y) plane, so the global coords are sequential
            int nzu = (nx * 2 + 1) * (ny + 1) + (nx + 1) * ny;
            var ZU = Enumerable.Range(0, nzu).ToArray();

#if DEBUG
            Debug.WriteLine("ZU: ");
            Debug.WriteLine(string.Join("; ", ZU));
#endif

            return ZU;
        }

        private int[,] CalculateZP(int nx, int ny, int nz)
        {
            // force is applied straight down at each vertext on the top plane
            int nep = (nx * 2 + 1) * (ny + 1) + (nx + 1) * ny;      // count of vertices on the top plane
            var ZP = new int[nep, 3];

            int totalVertexCount = ((nx * 2 + 1) + (nx + 1) * ny) * (nz + 1) + (nx + 1) * (ny + 1) * nz;
            int startPoint = totalVertexCount - nep;
            for (int i = 0; i < nep; i++)
            {
                ZP[i, 0] = startPoint + i;
                ZP[i, 1] = 5;                                       // 5 is top most (6th) plane of 1x1 cube
                ZP[i, 2] = 1;                                       // force applied
            }

#if DEBUG
            Debug.WriteLine("ZP: ");
            for (int i = 0; i < nep; i++)
            {
                Debug.WriteLine($"{ZP[i, 0]}  {ZP[i, 1]}  {ZP[i, 2]}");
            }
#endif

            return ZP;
        }
        
        // functions for dimensions for i < 8
        private Func<Point3D, Point3D, double>[] Dphis1 = new Func<Point3D, Point3D, double>[3]
        {
            (Point3D p, Point3D pi) => pi.X * (1 + p.Y * pi.Y) * (1 + p.Z * pi.Z) * (2 * p.X * pi.X + p.Y * pi.Y + p.Z * pi.Z - 1) / 8.0,
            (Point3D p, Point3D pi) => pi.Y * (1 + p.X * pi.X) * (1 + p.Z * pi.Z) * (2 * p.Y * pi.Y + p.X * pi.X + p.Z * pi.Z - 1) / 8.0,
            (Point3D p, Point3D pi) => pi.Z * (1 + p.X * pi.X) * (1 + p.Y * pi.Y) * (2 * p.Z * pi.Z + p.X * pi.X + p.Y * pi.Y - 1) / 8.0
        };

        // functions for dimensions for i >= 8
        private Func<Point3D, Point3D, double>[] Dphis2 = new Func<Point3D, Point3D, double>[3]
        {
            (Point3D p, Point3D pi) => (1 + p.Y * pi.Y) * (1 + p.Z * pi.Z) * (pi.Y * pi.Z * pi.Z + pi.X * pi.X * (pi.Y * pi.Y * p.Z + p.Y * pi.Z * pi.Z) + pi.X * (2 * pi.Y * p.X * pi.Z * pi.Z - 1)) / -4.0,
            (Point3D p, Point3D pi) => (1 + p.X * pi.X) * (1 + p.Z * pi.Z) * (pi.X * pi.Y * pi.Y * pi.Y * p.Z + pi.X * pi.Z * pi.Z + p.X * pi.Y * pi.Y * pi.Z * pi.Z + pi.Y * (2 * pi.X * p.Y * pi.Z * pi.Z - 1)) / -4.0,
            (Point3D p, Point3D pi) => (1 + p.X * pi.X) * (1 + p.Y * pi.Y) * (pi.X * p.Y * pi.Z * pi.Z * pi.Z + pi.X * pi.Y * pi.Y * (1 + 2 * p.Z * pi.Z) + pi.Z * (p.X * pi.Y * pi.Z * pi.Z - 1)) / -4.0
        };

        private double[,,] GenerateDFIABG()
        {
            var gaussPoints = Constants.GaussianAllCubePoints;

            // 20 points for alpha_i, beta_i and gamma_i
            var points = Constants.GaussianStandardCubePoints;
            
            // calculate DFIABG itself
            var result = new double[27,3,20];
            for (int cg = 0; cg < 27; ++cg)             // gauss points
            {
                var p = gaussPoints[cg];
                for (int d = 0; d < 3; ++d)             // dimension
                {
                    for (int i = 0; i < 20; ++i)        // functions
                    {
                        var funcArray = i < 8 ? Dphis1 : Dphis2;
                        result[cg, d, i] = funcArray[d](p, points[i]);
                    }
                }
            }

#if DEBUG
            Debug.WriteLine("DFIABG: ");
            for (int cg = 0; cg < 27; ++cg)
            {
                for (int d = 0; d < 3; ++d)
                {
                    for (int i = 0; i < 20; ++i)
                    {
                        Debug.WriteLine(result[cg, d, i]);
                    }
                }
            }
#endif

            return result;
        }
        
        private (double[,], double[]) ProcessElements(int nx, int ny, int nz, Point3D[] AKT, int[,] NT, double[,,] DFIABG, double[,,] DPSITE, int[,] ZP)
        {
            int npq = AKT.Length;
            var MG = new double[3 * npq, 3 * npq];
            var F = new double[3 * npq];

            int ce = nx * ny * nz;
            for (int i = 0; i < ce; ++i)
            {
#if DEBUG
                Debug.WriteLine($"Processing element #{i}");
#endif

                var DFIXYZ = CalculateDFIXYZ(i, AKT, NT);
                var DXYZABG = CalculateDXYZABG(i, AKT, NT, DFIABG);
                var DJ = CalculateDJ(DXYZABG);
                var MGE = CalculateMGE(DFIXYZ, DJ);
                var FE = CalculateFE(i, ce, DPSITE, ZP);

                UpdateMGF(MG, F, MGE, FE, NT, i);
            }

#if DEBUG
            Debug.WriteLine($"MG: ");
            for (int i = 0; i < 3 * npq; ++i)
            {
                for (int j = 0; j < 3 * npq; ++j)
                {
                    Debug.Write(MG[i, j] + " ");
                }

                Debug.WriteLine(string.Empty);
            }

            Debug.WriteLine($"F: ");
            for (int i = 0; i < 3 * npq; ++i)
            {
                Debug.WriteLine(F[i]);
            }
#endif

            return (MG, F);
        }

        private void UpdateMGF(double[,] MG, double[] F, double[,] MGE, double[] FE, int[,] NT, int feIndex)
        {
            for (int i = 0; i < 60; ++i)
            {
                int di = i / 20;                    // dimension (x,y or z ) for row
                int gi = NT[i % 20, feIndex];       // global index for row
                for (int j = 0; j < 60; ++j)
                {
                    int dj = j / 20;                // dimension (x,y or z ) for col
                    int gj = NT[j % 20, feIndex];   // global index for col

                    MG[3 * gi + di, 3 * gj + dj] = MGE[i, j];
                }

                F[3 * gi + di] = FE[i];
            }
        }

        private double[,,] CalculateDFIXYZ(int feIndex, Point3D[] AKT, int[,] NT)
        {
            var gaussPoints = Constants.GaussianAllCubePoints;

            // calculate DFIXYZ itself
            var result = new double[27, 20, 3];
            for (int cg = 0; cg < 27; ++cg)             // gauss points
            {
                var p = gaussPoints[cg];
                for (int d = 0; d < 3; ++d)             // dimension
                {
                    for (int i = 0; i < 20; ++i)        // functions
                    {
                        var pi = AKT[NT[i, feIndex]];   // use NT to lookup global point using feIndex and i
                        var funcArray = i < 8 ? Dphis1 : Dphis2;
                        result[cg, i, d] = funcArray[d](p, pi);
                    }
                }
            }

#if DEBUG
            Debug.WriteLine($"DFIXYZ (feIndex = {feIndex}): ");
            for (int cg = 0; cg < 27; ++cg)
            {
                for (int d = 0; d < 3; ++d)
                {
                    for (int i = 0; i < 20; ++i)
                    {
                        Debug.WriteLine(result[cg, i, d]);
                    }
                }
            }
#endif

            return result;
        }


        private double[,,] CalculateDXYZABG(int feIndex, Point3D[] AKT, int[,] NT, double[,,] DFIABG)
        {
            var gaussPoints = Constants.GaussianAllCubePoints;

            Func<Point3D, double>[] globalValuesGetters =
            {
                p => p.X,
                p => p.Y,
                p => p.Z
            };

            // calculate DXYZABG itself
            var result = new double[3, 3, 27];
            for (int gd = 0; gd < 3; ++gd)                 // global coord
            {
                var p = gaussPoints[gd];
                for (int ld = 0; ld < 3; ++ld)             // local coord
                {
                    for (int cg = 0; cg < 27; ++cg)        // gauss points
                    {
                        var getter = globalValuesGetters[gd];
                        double localSum = 0.0;
                        for (int i = 0; i < 20; ++i)       // functions
                        {
                            var pi = AKT[NT[i, feIndex]];  // global point
                            double piValue = getter(pi);   // value of global coordinate
                            double dphi = DFIABG[cg, ld, i];
                            localSum += piValue * dphi;
                        }

                        result[gd, ld, cg] = localSum;
                    }
                }
            }

#if DEBUG
            Debug.WriteLine($"DXYZABG (feIndex = {feIndex}): ");
            for (int gd = 0; gd < 3; ++gd)
            {
                for (int ld = 0; ld < 3; ++ld)
                {
                    for (int cg = 0; cg < 27; ++cg)
                    {
                        Debug.WriteLine(result[gd, ld, cg]);
                    }
                }
            }
#endif

            return result;
        }

        private double Matrix3by3Det(double a11, double a12, double a13
                                   , double a21, double a22, double a23
                                   , double a31, double a32, double a33)
        {
            return a11 * a22 * a33 + a12 * a23 * a31 + a13 * a21 * a32
                 - a13 * a22 * a31 - a23 * a32 * a11 - a33 * a12 * a21;
        }

        private double[] CalculateDJ(double[,,] DXYZABG)
        {
            var result = new double[27];            // length == count of gaussian points
            for (int cg = 0; cg < 27; ++cg)         // gaussian points
            {
                result[cg] = Matrix3by3Det(DXYZABG[0, 0, cg], DXYZABG[1, 0, cg], DXYZABG[2, 0, cg]
                                         , DXYZABG[0, 1, cg], DXYZABG[1, 1, cg], DXYZABG[2, 1, cg]
                                         , DXYZABG[0, 2, cg], DXYZABG[1, 2, cg], DXYZABG[2, 2, cg]);
            }
            
#if DEBUG
            Debug.WriteLine($"DJ: ");
            for (int cg = 0; cg < 27; ++cg)
            {
                Debug.WriteLine(result[cg]);
            }
#endif

            return result;
        }

        private double[,] CalculateMGE(double[,,] DFIXYZ, double[] DJ)
        {
            var Cs = Constants.Cs;
            Func<int, int, int, double>[,] As =
            {
                {
                    (int i, int j, int cg) => l * (1.0 - v * DFIXYZ[cg, i, 0] * DFIXYZ[cg, j, 0]) + m * (DFIXYZ[cg, i, 1] * DFIXYZ[cg, j, 1] + DFIXYZ[cg, i, 2] * DFIXYZ[cg, j, 2]),
                    (int i, int j, int cg) => l * v * DFIXYZ[cg, i, 0] * DFIXYZ[cg, j, 1] + m * DFIXYZ[cg, i, 1] * DFIXYZ[cg, j, 0],
                    (int i, int j, int cg) => l * v * DFIXYZ[cg, i, 0] * DFIXYZ[cg, j, 2] + m * DFIXYZ[cg, i, 2] * DFIXYZ[cg, j, 0]
                },
                {
                    null,
                    (int i, int j, int cg) => l * (1.0 - v * DFIXYZ[cg, i, 1] * DFIXYZ[cg, j, 1]) + m * (DFIXYZ[cg, i, 0] * DFIXYZ[cg, j, 0] + DFIXYZ[cg, i, 2] * DFIXYZ[cg, j, 2]),
                    (int i, int j, int cg) => l * v * DFIXYZ[cg, i, 1] * DFIXYZ[cg, j, 2] + m * DFIXYZ[cg, i, 2] * DFIXYZ[cg, j, 1]
                },
                {
                    null,
                    null,
                    (int i, int j, int cg) => l * (1.0 - v * DFIXYZ[cg, i, 2] * DFIXYZ[cg, j, 2]) + m * (DFIXYZ[cg, i, 0] * DFIXYZ[cg, j, 0] + DFIXYZ[cg, i, 1] * DFIXYZ[cg, j, 1])
                }
            };
            As[1, 0] = As[0, 1];
            As[2, 0] = As[0, 2];
            As[2, 1] = As[1, 2];

            var result = new double[60, 60];
            for (int ai = 0; ai < 3; ++ai)
            {
                for (int aj = 0; aj < 3; ++aj)
                {
                    var a = As[ai, aj];
                    for (int i = 0; i < 20; ++i)
                    {
                        for (int j = 0; j < 20; ++j)
                        {
                            int cg = 0;
                            double res = 0.0;
                            for (int m = 0; m < 3; ++m)
                            {
                                for (int n = 0; n < 3; ++n)
                                {
                                    for (int k = 0; k < 3; ++k)
                                    {
                                        res += Cs[m] * Cs[n] * Cs[k] * a(i, j, cg) * DJ[cg];
                                        ++cg;
                                    }
                                }
                            }

                            result[ai * 20 + i, aj * 20 + j] = res;
                        }
                    }
                }
            }

#if DEBUG
            Debug.WriteLine($"MGE: ");
            for (int i = 0; i < 60; ++i)
            {
                for (int j = 0; j < 60; ++j)
                {
                    Debug.WriteLine(result[i, j]);
                }
            }
#endif

            return result;
        }
        
        // TODO: calculate FE
        private double[] CalculateFE(int feIndex, int feCount, double[,,] DPSITE, int[,] ZP)
        {
            var result = new double[60];

#if DEBUG
            Debug.WriteLine($"FE: ");
            for (int i = 0; i < 60; ++i)
            {
                Debug.WriteLine(result[i]);
            }
#endif

            return result;
        }
        
        private void FixMG(double[,] MG, int[] ZU)
        {
            for (int d = 0; d < 3; ++d)
            {
                foreach (var i in ZU)
                {
                    int index = 3 * i + d;
                    MG[index, index] = LargeCoefficient;
                }
            }
        }

        private double[] GaussianElimination(double[,] A, double[] b)
        {
            int n = b.Length;
            var x = new double[n];

            // append b to the right of A
            var AB = new double[n, n + 1];
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    AB[i, j] = A[i, j];
                }

                AB[i, n] = b[i];
            }

            for (int i = 0; i < n; i++)
            {
                // Search for maximum in this column
                double maxEl = Math.Abs(AB[i, i]);
                int maxRow = i;
                for (int k = i + 1; k < n; k++)
                {
                    if (Math.Abs(AB[k, i]) > maxEl)
                    {
                        maxEl = Math.Abs(AB[k, i]);
                        maxRow = k;
                    }
                }

                // Swap maximum row with current row (column by column)
                for (int k = i; k < n + 1; k++)
                {
                    double tmp = AB[maxRow, k];
                    AB[maxRow, k] = AB[i, k];
                    AB[i, k] = tmp;
                }

                // Make all rows below this one 0 in current column
                for (int k = i + 1; k < n; k++)
                {
                    double c = -AB[k, i] / AB[i, i];
                    for (int j = i; j < n + 1; j++)
                    {
                        if (i == j)
                        {
                            AB[k, j] = 0;
                        }
                        else
                        {
                            AB[k, j] += c * AB[i, j];
                        }
                    }
                }
            }

            // Solve equation Ax=b for an upper triangular matrix A
            for (int i = n - 1; i >= 0; i--)
            {
                x[i] = AB[i, n] / AB[i, i];
                for (int k = i - 1; k >= 0; k--)
                {
                    AB[k, n] -= AB[k, i] * x[i];
                }
            }

            return x;
        }

        private Point3D[] AddTranslation(Point3D[] AKT, double[] U)
        {
            var result = new Point3D[AKT.Length];
            for (int i = 0; i < AKT.Length; ++i)
            {
                result[i] = new Point3D()
                {
                    X = AKT[i].X + U[i * 3],
                    Y = AKT[i].Y + U[i * 3 + 1],
                    Z = AKT[i].Z + U[i * 3 + 2]
                };
            }

            return result;
        }

        private void GenerateAndAddShapeToCollection(Point3D[] points, int[] pointIndexes, Brush vertColor, Visual3DCollection collection)
        {
            // add all points (vertices)
            foreach (var point in points)
            {
                var vert = CreateVertex(point, vertColor);
                collection.Add(vert);
            }

            // TODO: add lines
            // add all lines (edges)
            //for (int i = 0; i < pointIndexes.Length; i += 2)
            //{
            //    var line = CreateLine(points[pointIndexes[i]], points[pointIndexes[i + 1]]);
            //    collection.Add(line);
            //}
        }

        private PipeVisual3D CreateLine(Point3D start, Point3D end)
        {
            var line = new PipeVisual3D();

            line.Point1 = start;
            line.Point2 = end;
            line.InnerDiameter = LineInnerDiameter;
            line.Diameter = LineDiameter;
            line.Fill = LineFill;

            return line;
        }

        private SphereVisual3D CreateVertex(Point3D center, Brush fill)
        {
            var vert = new SphereVisual3D();

            vert.Center = center;
            vert.Radius = VertRadius;
            vert.Fill = fill;

            return vert;
        }
    }
}
