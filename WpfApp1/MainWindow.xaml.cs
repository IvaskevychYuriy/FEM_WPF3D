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
        public MainWindow()
        {
            InitializeComponent();
            Render();
        }

        private void Render()
        {
            var points = new List<Point3D>()
            {
                new Point3D(-1, 0, 0), new Point3D(0, 0, 0), new Point3D(0, -1, 0), new Point3D(-1, -1, 0)
            };

            var edgesInfo = new List<int>()
            {
                0, 1, 1, 2, 2, 3, 3, 0
            };

            GenerateAndAddShapeToCollection(points, edgesInfo, Brushes.Blue, MainViewPort.Children);
        }

        private void GenerateAndAddShapeToCollection(List<Point3D> points, List<int> pointIndexes, Brush vertColor, Visual3DCollection collection)
        {
            // add all points (vertices)
            foreach (var point in points)
            {
                var vert = CreateVertex(point, vertColor);
                collection.Add(vert);
            }

            // add all lines (edges)
            for (int i = 0; i < pointIndexes.Count; i += 2)
            {
                var line = CreateLine(points[pointIndexes[i]], points[pointIndexes[i+1]]);
                collection.Add(line);
            }
        }

        private PipeVisual3D CreateLine(Point3D start, Point3D end)
        {
            var line = new PipeVisual3D();

            line.Point1 = start;
            line.Point2 = end;
            line.InnerDiameter = 0;
            line.Diameter = 0.05;
            line.Fill = Brushes.LightGray;

            return line;
        }

        private SphereVisual3D CreateVertex(Point3D center, Brush fill)
        {
            var vert = new SphereVisual3D();

            vert.Center = center;
            vert.Radius = 0.1;
            vert.Fill = fill;

            return vert;
        }
    }
}
