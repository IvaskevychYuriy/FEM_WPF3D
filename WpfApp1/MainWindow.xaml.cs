﻿using HelixToolkit.Wpf;
using MathLib;
using System;
using System.Collections.Generic;
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
        private const double v = 0.3;
        private const double l = E / ((1 + v) * (1 - 2 * v));
        private const double m = E / (2 * (1 + v));
        private const double Pressure = -0.1;

        private readonly double LargeCoefficient = Math.Pow(10, 20);

        public MainWindow()
        {
            InitializeComponent();

            const int nx = 1;
            const int ny = 1;
            const int nz = 1;
                           
            const int ax = 1;
            const int ay = 1;
            const int az = 1;

            // initial points array
            var AKT = GenerateAKT(nx, ny, nz, ax, ay, az);
            Render(AKT);
            
            var NT = GeneratorNP.Generate(globalToAKT, nx, ny, nz);

            // fixed points global coords (ZU) and force info (ZP)
            var ZU = CalculateZU(nx, ny, nz);
            var ZP = CalculateZP(nx, ny, nz, AKT.Length);

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
            WriteLine("Initial points (AKT): ");
            foreach (var p in points)
            {
                WriteLine(p);
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
            WriteLine("ZU: ");
            WriteLine(string.Join("; ", ZU));
#endif

            return ZU;
        }

        private double[,] CalculateZP(int nx, int ny, int nz, int totalVertexCount)
        {
            // force is applied straight down at each vertext on the top plane
            int nep = (nx * 2 + 1) * (ny + 1) + (nx + 1) * ny;      // count of vertices on the top plane
            var ZP = new double[nep, 3];
            
            int startPoint = totalVertexCount - nep;
            for (int i = 0; i < nep; i++)
            {
                ZP[i, 0] = startPoint + i;
                ZP[i, 1] = 5;                                       // 5 is top most (6th) plane of 1x1 cube
                ZP[i, 2] = Pressure;
            }

#if DEBUG
            WriteLine("ZP: ");
            for (int i = 0; i < nep; i++)
            {
                WriteLine($"{ZP[i, 0]}  {ZP[i, 1]}  {ZP[i, 2]}");
            }
#endif

            return ZP;
        }

        // functions for dimensions for i < 8
        private Func<Point3D, Point3D, double>[] Dphis1 = new Func<Point3D, Point3D, double>[3]
        {
            (Point3D p, Point3D pi) => 0.125 * pi.X * (1 + p.Y * pi.Y) * (1 + p.Z * pi.Z) * (2 * p.X * pi.X + p.Y * pi.Y + p.Z * pi.Z - 1),
            (Point3D p, Point3D pi) => 0.125 * pi.Y * (1 + p.X * pi.X) * (1 + p.Z * pi.Z) * (2 * p.Y * pi.Y + p.X * pi.X + p.Z * pi.Z - 1),
            (Point3D p, Point3D pi) => 0.125 * pi.Z * (1 + p.X * pi.X) * (1 + p.Y * pi.Y) * (2 * p.Z * pi.Z + p.X * pi.X + p.Y * pi.Y - 1)
        };

        // functions for dimensions for i >= 8
        private Func<Point3D, Point3D, double>[] Dphis2 = new Func<Point3D, Point3D, double>[3]
        {
            (Point3D p, Point3D pi) =>
            0.25 *
                (1 + pi.Y * p.Y) *
                (1 + pi.Z * p.Z) *
                (pi.X * (1 - Math.Pow((p.X * pi.Y * pi.Z), 2) - Math.Pow((p.Y * pi.X * pi.Z), 2) - Math.Pow((p.Z * pi.X * pi.Y), 2)) -
                  2 * p.X * Math.Pow((pi.Y * pi.Z), 2) * (1 + p.X * pi.X))
            
            ,
            (Point3D p, Point3D pi) =>
            0.25 *
                (1 + pi.X * p.X) *
                (1 + pi.Z * p.Z) *
                (pi.Y * (1 - Math.Pow((p.X * pi.Y * pi.Z), 2) - Math.Pow((p.Y * pi.X * pi.Z), 2) - Math.Pow((p.Z * pi.X * pi.Y), 2)) -
                  2 * p.Y * Math.Pow((pi.X * pi.Z), 2) * (1 + p.Y * pi.Y))
            ,
            (Point3D p, Point3D pi) =>
            0.25 *
                (1 + pi.X * p.X) *
                (1 + pi.Y * p.Y) *
                (pi.Z * (1 - Math.Pow((p.X * pi.Y * pi.Z), 2) - Math.Pow((p.Y * pi.X * pi.Z), 2) - Math.Pow((p.Z * pi.X * pi.Y), 2)) -
                  2 * p.Z * Math.Pow((pi.X * pi.Y), 2) * (1 + p.Z * pi.Z))
        };

        private double[,,] GenerateDFIABG()
        {
            var gaussPoints = Constants.GaussianAllCubePoints;

            // 20 points for alpha_i, beta_i and gamma_i
            var points = Constants.GaussianStandardCubePoints;

            // calculate DFIABG itself
            var result = new double[27, 3, 20];
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
            WriteLine("DFIABG: ");
            for (int cg = 0; cg < 27; ++cg)
            {
                for (int d = 0; d < 3; ++d)
                {
                    for (int i = 0; i < 20; ++i)
                    {
                        WriteLine(result[cg, d, i]);
                    }
                }
            }
#endif

            return result;
        }

        private (double[,], double[]) ProcessElements(int nx, int ny, int nz, Point3D[] AKT, int[,] NT, double[,,] DFIABG, double[,,] DPSITE, double[,] ZP)
        {
            int npq = AKT.Length;
            var MG = new double[3 * npq, 3 * npq];
            var F = new double[3 * npq];

            int ce = nx * ny * nz;
            for (int i = 0; i < ce; ++i)
            {
#if DEBUG
                WriteLine($"Processing element #{i}");
#endif

                var DXYZABG = CalculateDXYZABG(i, AKT, NT, DFIABG);
                var DFIXYZ = CalculateDFIXYZ(i, AKT, NT, DFIABG, DXYZABG);
                var DJ = CalculateDJ(DXYZABG);
                var MGE = CalculateMGE(DFIXYZ, DJ);
                var FE = CalculateFE(i, ce, nx * ny, DPSITE, ZP, NT, AKT);

                UpdateMGF(MG, F, MGE, FE, NT, i);
            }

#if DEBUG
            WriteLine($"MG: ");
            for (int i = 0; i < 3 * npq; ++i)
            {
                for (int j = 0; j < 3 * npq; ++j)
                {
                    Write(MG[i, j] + " ");
                }

                WriteLine(string.Empty);
            }

            WriteLine($"F: ");
            for (int i = 0; i < 3 * npq; ++i)
            {
                WriteLine(F[i]);
            }
#endif

            return (MG, F);
        }

        private void UpdateMGF(double[,] MG, double[] F, double[,] MGE, double[] FE, int[,] NT, int feIndex)
        {
            var pressedLocalPoints = Constants.PressedLocalPoints;

            // update MG
            for (int i = 0; i < 60; ++i)
            {
                int di = i / 20;                    // dimension (x,y or z) for row
                int gi = NT[i % 20, feIndex];       // global index for row
                for (int j = 0; j < 60; ++j)
                {
                    int dj = j / 20;                // dimension (x,y or z) for col
                    int gj = NT[j % 20, feIndex];   // global index for col

                    MG[3 * gi + di, 3 * gj + dj] += MGE[i, j];
                }
            }

            // update F
            int startIndex = 60 - 24 + 2;           // third (for Z-coord) starting from last 24 in all FE (60)
            for (int i = 0; i < pressedLocalPoints.Length; ++i)
            {
                F[NT[pressedLocalPoints[i], feIndex] * 3 + 2] += FE[startIndex]; // z-coords only
                startIndex += 3;
            }
        }

        private double[,,] CalculateDFIXYZ(int feIndex, Point3D[] AKT, int[,] NT, double[,,] DFIABG, double[,,] DXYZABG)
        {
            var gaussPoints = Constants.GaussianAllCubePoints;

            // calculate DFIXYZ itself
            var result = new double[27, 20, 3];
            for (int cg = 0; cg < 27; ++cg)             // gauss points
            {
                var p = gaussPoints[cg];          
                {
                    for (int i = 0; i < 20; ++i)        // functions
                    {
                        var b = new double[] { DFIABG[cg, 0, i], DFIABG[cg, 1, i], DFIABG[cg, 2, i] };
                        var res = GaussianElimination(CalculateD(DXYZABG, cg), b);
                        for (int d = 0; d < 3; ++d)     // dimension
                        {
                            result[cg, i, d] = res[d];
                        }
                    }
                }
            }

#if DEBUG
            WriteLine($"DFIXYZ (feIndex = {feIndex}): ");
            for (int cg = 0; cg < 27; ++cg)
            {
                for (int d = 0; d < 3; ++d)
                {
                    for (int i = 0; i < 20; ++i)
                    {
                        WriteLine(result[cg, i, d]);
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
                var getter = globalValuesGetters[gd];
                for (int ld = 0; ld < 3; ++ld)             // local coord
                {
                    for (int cg = 0; cg < 27; ++cg)        // gauss points
                    {
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
            WriteLine($"DXYZABG (feIndex = {feIndex}): ");
            for (int gd = 0; gd < 3; ++gd)
            {
                for (int ld = 0; ld < 3; ++ld)
                {
                    for (int cg = 0; cg < 27; ++cg)
                    {
                        WriteLine(result[gd, ld, cg]);
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
            WriteLine($"DJ: ");
            for (int cg = 0; cg < 27; ++cg)
            {
                WriteLine(result[cg]);
            }
#endif

            return result;
        }

        private double[,] CalculateD(double[,,] DXYZABG, int cg)
        {
            var result = new double[,]
            {
                { DXYZABG[0, 0, cg], DXYZABG[1, 0, cg], DXYZABG[2, 0, cg] },
                { DXYZABG[0, 1, cg], DXYZABG[1, 1, cg], DXYZABG[2, 1, cg] },
                { DXYZABG[0, 2, cg], DXYZABG[1, 2, cg], DXYZABG[2, 2, cg] }
            };

#if DEBUG
            WriteLine($"D (cg={cg}): ");
            for (int i = 0; i < 3; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    Write(result[i,j] + " ");
                }

                WriteLine("");
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
                    (int i, int j, int cg) => (l * (1.0 - v) * DFIXYZ[cg, i, 0] * DFIXYZ[cg, j, 0]) + m * (DFIXYZ[cg, i, 1] * DFIXYZ[cg, j, 1] + DFIXYZ[cg, i, 2] * DFIXYZ[cg, j, 2]),
                    (int i, int j, int cg) => l * v * DFIXYZ[cg, i, 0] * DFIXYZ[cg, j, 1] + m * DFIXYZ[cg, i, 1] * DFIXYZ[cg, j, 0],
                    (int i, int j, int cg) => l * v * DFIXYZ[cg, i, 0] * DFIXYZ[cg, j, 2] + m * DFIXYZ[cg, i, 2] * DFIXYZ[cg, j, 0]
                },
                {
                    null,
                    (int i, int j, int cg) => (l * (1.0 - v) * DFIXYZ[cg, i, 1] * DFIXYZ[cg, j, 1]) + m * (DFIXYZ[cg, i, 0] * DFIXYZ[cg, j, 0] + DFIXYZ[cg, i, 2] * DFIXYZ[cg, j, 2]),
                    (int i, int j, int cg) => l * v * DFIXYZ[cg, i, 1] * DFIXYZ[cg, j, 2] + m * DFIXYZ[cg, i, 2] * DFIXYZ[cg, j, 1]
                },
                {
                    null,
                    null,
                    (int i, int j, int cg) => (l * (1.0 - v) * DFIXYZ[cg, i, 2] * DFIXYZ[cg, j, 2]) + m * (DFIXYZ[cg, i, 0] * DFIXYZ[cg, j, 0] + DFIXYZ[cg, i, 1] * DFIXYZ[cg, j, 1])
                }
            };

            var result = new double[60, 60];
            for (int ai = 0; ai < 3; ++ai)
            {
                for (int aj = 0; aj < 3; ++aj)
                {
                    var a = As[ai, aj];
                    if (a == null)
                    {
                        continue;
                    }

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
                                        res += Cs[m] * Cs[n] * Cs[k] * a(i, j, cg) * Math.Abs(DJ[cg]);
                                        ++cg;
                                    }
                                }
                            }

                            result[ai * 20 + i, aj * 20 + j] = res;
                        }
                    }
                }
            }

            for (int i = 0; i < 60; ++i)
            {
                for (int j = 0; j < i; ++j)
                {
                    result[i, j] = result[j, i];
                }
            }

#if DEBUG
            bool symmetric = true;
            bool diag = true;
            for (int i = 0; i < 60; ++i)
            {
                if (result[i,i] < 0)
                {
                    diag = false;
                }
                for (int j = 0; j < i; ++j)
                {
                    if (result[i, j] != result[j, i])
                    {
                        symmetric = false;
                    }
                }
            }
            
            WriteLine($"MGE ({symmetric} {diag}): ");
            for (int i = 0; i < 60; ++i)
            {
                for (int j = 0; j < 60; ++j)
                {
                    WriteLine(result[i, j]);
                }
            }
#endif

            return result;
        }
        
        private double[] CalculateFE(int feIndex, int feCount, int feCountUnderPressure, double[,,] DPSITE, double[,] ZP, int[,] NT, Point3D[] AKT)
        {
            var result = new double[60];
            if (feIndex < feCount - feCountUnderPressure)
            {
                // current element not under force
                return result;
            }

            var nodes = GenerateDPSITE.GetNiTi();
            var dPsiFunctors = new Func<int, double>[]
            {
                (i) => GenerateDPSITE.DN(i, nodes[i].Item1, nodes[i].Item2),
                (i) => GenerateDPSITE.DT(i, nodes[i].Item1, nodes[i].Item2),
            };

            var globalCoordGetters = new Func<Point3D, double>[]
            {
                p => p.X,
                p => p.Y,
                p => p.Z
            };

            var pressedLocalPoints = Constants.PressedLocalPoints;
            var gaussianNodes = GenerateDPSITE.GetGaussNode();
            var derivatives = new double[3, 2, 9];                  // dxyz / dnt
            for (int d1 = 0; d1 < 3; ++d1)                          // 1st dimention (x, y, z)
            {
                for (int d2 = 0; d2 < 2; ++d2)                      // 2nd dimention (n, t)
                {
                    for (int gp = 0; gp < 9; gp++)
                    {
                        double localSum = 0.0;
                        for (int i = 0; i < 8; ++i)
                        {
                            var globalPoint = AKT[NT[pressedLocalPoints[i], feIndex]];
                            double coordValue = globalCoordGetters[d1](globalPoint);
                            localSum += coordValue * DPSITE[gp, d2, i];
                        }

                        derivatives[d1, d2, gp] = localSum;
                    }
                }
            }

            var zpMap = new Dictionary<double, int>();
            for (int i = 0; i < ZP.GetLength(0); i++)
            {
                zpMap.Add(ZP[i, 0], i);                    // global point index - index in ZP
            }

            int startIndex = 60 - 24 + 2;                  // third (for Z-coord) starting from last 24 in all FE (60)
            var Cs = Constants.Cs;
            for (int i = 0; i < 8; ++i)
            {
                double localSum = 0.0;
                for (int m = 0; m < 3; ++m)
                {
                    for (int n = 0; n < 3; ++n)
                    {
                        int cg = m * 3 + n;
                        localSum += Cs[m] * Cs[n] * Pressure * (derivatives[0, 0, cg] * derivatives[1, 1, cg] - derivatives[1, 0, cg] * derivatives[0, 1, cg]) * GenerateDPSITE.PHIs[i](gaussianNodes[cg, 0], gaussianNodes[cg, 1]);
                    }
                }

                result[startIndex] = localSum;
                startIndex += 3;
            }

#if DEBUG
            WriteLine($"FE: ");
            for (int i = 0; i < 60; ++i)
            {
                WriteLine(result[i]);
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

        private static string Id { get; } = DateTime.Now.ToString().Replace(" ", "_").Replace(":", "-");

        private void WriteLine<T>(T str)
        {
            using (var s = new System.IO.StreamWriter(System.IO.Path.Combine("../../", "test-" + Id + ".txt"), true))
            {
                s.WriteLine(str);
                s.Flush();
            }
        }

        private void Write(string str)
        {
            using (var s = new System.IO.StreamWriter(System.IO.Path.Combine("../../", "test-" + Id + ".txt"), true))
            {
                s.Write(str);
                s.Flush();
            }
        }
    }
}
