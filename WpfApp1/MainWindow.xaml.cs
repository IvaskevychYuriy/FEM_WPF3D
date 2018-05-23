using HelixToolkit.Wpf;
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

        public MainWindow()
        {
            InitializeComponent();
         
            const int nx = 2;
            const int ny = 1;
            const int nz = 2;

            // initial points array
            var AKT = GenerateAKT(nx, ny, nz);
            Render(AKT);

            // TODO: calculate NT
            // var NT = ...

            // fixed points global coords (ZU) and force info (ZP)
            var ZU = CalculateZU(nx, ny, nz);
            var ZP = CalculateZP(nx, ny, nz);

            var DFIABG = GenerateDFIABG();
        }

        private void Render(Point3D[] points)
        {
            // TODO: add edges
            GenerateAndAddShapeToCollection(points, new int[] { }, Brushes.Blue, MainViewPort.Children);
        }

        private Point3D[] GenerateAKT(int nx, int ny, int nz)
        {
            int cx = nx + 1;
            int cy = ny + 1;
            int cz = nz + 1;
            int cxyz = cx * cy * cz;                          // 1x1 boxes vertex count
            int ce = (cx * nz + nx * cz) * cy + ny * cx * cz; // edges count

            var points = new Point3D[cxyz + ce];

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
                            points[i++] = new Point3D(ix / 2.0, iy / 2.0, iz / 2.0);
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

        private double[,,] GenerateDFIABG()
        {
            // generate 27 points in -1 to 1 standard cube 
            var gaussPoints = new Point3D[27];
            int c = 0;
            for (int cx = 0; cx < 3; ++cx)
            {
                for (int cy = 0; cy < 3; ++cy)
                {
                    for (int cz = 0; cz < 3; ++cz)
                    {
                        gaussPoints[c++] = new Point3D(-1.0 + cx * 0.5, -1.0 + cy * 0.5, -1.0 + cz * 0.5);
                    }
                }
            }

            // 20 points for alpha_i, beta_i and gamma_i
            var points = new Point3D[20]
            {
                // first 8
                new Point3D(-1, -1, -1),
                new Point3D(1, -1, -1),
                new Point3D(1, 1, -1),
                new Point3D(-1, 1, -1),
                new Point3D(-1, -1, 1),
                new Point3D(1, -1, 1),
                new Point3D(1, 1, 1),
                new Point3D(-1, 1, 1),

                // rest 12
                new Point3D(0, -1, -1),
                new Point3D(1, 0, -1),
                new Point3D(0, 1, -1),
                new Point3D(-1, 0, -1),
                new Point3D(-1, -1, 0),
                new Point3D(1, -1, 0),
                new Point3D(1, 1, 0),
                new Point3D(-1, 1, 0),
                new Point3D(0, -1, 1),
                new Point3D(1, 0, 1),
                new Point3D(0, 1, 1),
                new Point3D(-1, 0, 1)
            };

            // functions for dimensions for i < 8
            var dphis1 = new Func<Point3D, Point3D, double>[3]
            {
                (Point3D p, Point3D pi) => pi.X * (1 + p.Y * pi.Y) * (1 + p.Z * pi.Z) * (2 * p.X * pi.X + p.Y * pi.Y + p.Z * pi.Z - 1) / 8.0,
                (Point3D p, Point3D pi) => pi.Y * (1 + p.X * pi.X) * (1 + p.Z * pi.Z) * (2 * p.Y * pi.Y + p.X * pi.X + p.Z * pi.Z - 1) / 8.0,
                (Point3D p, Point3D pi) => pi.Z * (1 + p.X * pi.X) * (1 + p.Y * pi.Y) * (2 * p.Z * pi.Z + p.X * pi.X + p.Y * pi.Y - 1) / 8.0
            };

            // functions for dimensions for i >= 8
            var dphis2 = new Func<Point3D, Point3D, double>[3]
            {
                (Point3D p, Point3D pi) => (1 + p.Y * pi.Y) * (1 + p.Z * pi.Z) * (pi.Y * pi.Z * pi.Z + pi.X * pi.X * (pi.Y * pi.Y * p.Z + p.Y * pi.Z * pi.Z) + pi.X * (2 * pi.Y * p.X * pi.Z * pi.Z - 1)) / -4.0,
                (Point3D p, Point3D pi) => (1 + p.X * pi.X) * (1 + p.Z * pi.Z) * (pi.X * pi.Y * pi.Y * pi.Y * p.Z + pi.X * pi.Z * pi.Z + p.X * pi.Y * pi.Y * pi.Z * pi.Z + pi.Y * (2 * pi.X * p.Y * pi.Z * pi.Z - 1)) / -4.0,
                (Point3D p, Point3D pi) => (1 + p.X * pi.X) * (1 + p.Y * pi.Y) * (pi.X * p.Y * pi.Z * pi.Z * pi.Z + pi.X * pi.Y * pi.Y * (1 + 2 * p.Z * pi.Z) + pi.Z * (p.X * pi.Y * pi.Z * pi.Z - 1)) / -4.0
            };

            // calculate DFIABG itself
            var result = new double[27,3,20];
            for (int cg = 0; cg < 27; ++cg)             // gauss points
            {
                var p = gaussPoints[cg];
                for (int d = 0; d < 3; ++d)             // dimension
                {
                    for (int i = 0; i < 8; ++i)        // functions < 8
                    {
                        result[cg, d, i] = dphis1[d](p, points[i]);
                    }
                    for (int i = 8; i < 20; ++i)        // functions >= 8
                    {
                        result[cg, d, i] = dphis2[d](p, points[i]);
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
