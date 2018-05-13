using HelixToolkit.Wpf;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Media.Media3D;
using System.Windows.Navigation;
using System.Windows.Shapes;

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

        public MainWindow()
        {
            InitializeComponent();
            Render();
        }

        private void Render()
        {
            var data = GenerateInitialPoints(2, 2, 2);
            GenerateAndAddShapeToCollection(data.Vertexes, data.EdgesData, Brushes.Blue, MainViewPort.Children);
        }

        private GridData GenerateInitialPoints(int nx, int ny, int nz)
        {
            int cx = nx + 1;
            int cy = ny + 1;
            int cz = nz + 1;
            int cxyz = cx * cy * cz;                          // 1x1 boxes vertex count
            int ce = (cx * nz + nx * cz) * cy + ny * cx * cz; // edges count

            var points = new Point3D[cxyz + ce];
            var edgesData = new int[ce * 4];
            Func<int, int, int, int> calcIndex = (int x, int y, int z) => x + cx * (y + cy * z);

            int ic = 0;
            for (int ix = 0; ix < cx; ++ix)
            {
                for (int iy = 0; iy < cy; ++iy)
                {
                    for (int iz = 0; iz < cz; ++iz)
                    {
                        int i = calcIndex(ix, iy, iz);
                        points[i] = new Point3D(ix, iy, iz);

                        if (iz > 0)
                        {
                            points[cxyz + ic] = new Point3D(ix, iy, iz - 0.5);
                            int iprev = calcIndex(ix, iy, iz - 1);
                            FillEdgesData(cxyz, i, iprev, ref ic, edgesData);
                        }
                        if (iy > 0)
                        {
                            points[cxyz + ic] = new Point3D(ix, iy - 0.5, iz);
                            int iprev = calcIndex(ix, iy - 1, iz);
                            FillEdgesData(cxyz, i, iprev, ref ic, edgesData);
                        }
                        if (ix > 0)
                        {
                            points[cxyz + ic] = new Point3D(ix - 0.5, iy, iz);
                            int iprev = calcIndex(ix - 1, iy, iz);
                            FillEdgesData(cxyz, i, iprev, ref ic, edgesData);
                        }
                    }
                }
            }

            return new GridData()
            {
                Vertexes = points,
                EdgesData = edgesData
            };
        }
        
        private void FillEdgesData(int cxyz, int i, int iprev, ref int ic, int[] edgesData)
        {
            edgesData[ic * 4] = iprev;
            edgesData[ic * 4 + 1] = edgesData[ic * 4 + 2] = cxyz + ic;
            edgesData[ic * 4 + 3] = i;
            ++ic;
        }

        private void GenerateAndAddShapeToCollection(Point3D[] points, int[] pointIndexes, Brush vertColor, Visual3DCollection collection)
        {
            // add all points (vertices)
            foreach (var point in points)
            {
                var vert = CreateVertex(point, vertColor);
                collection.Add(vert);
            }

            // add all lines (edges)
            for (int i = 0; i < pointIndexes.Length; i += 2)
            {
                var line = CreateLine(points[pointIndexes[i]], points[pointIndexes[i + 1]]);
                collection.Add(line);
            }
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

    internal class GridData
    {
        public Point3D[] Vertexes { get; set; }
        public int[] EdgesData { get; set; }
    }
}
