using System.Windows.Media.Media3D;

namespace MathLib
{
    public static class Constants
    {
        public static Point3D[] GaussianAllCubePoints { get; }

        public static Point3D[] GaussianStandardCubePoints { get; } = new Point3D[20]
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

        static Constants()
        {
            GaussianAllCubePoints = GenerateGaussianPoints();
        }

        private static Point3D[] GenerateGaussianPoints()
        {
            // generate 27 points in -1 to 1 standard cube 
            var gaussianPoints = new Point3D[27];
            int c = 0;
            for (int cx = 0; cx < 3; ++cx)
            {
                for (int cy = 0; cy < 3; ++cy)
                {
                    for (int cz = 0; cz < 3; ++cz)
                    {
                        gaussianPoints[c++] = new Point3D(-1.0 + cx * 1.0, -1.0 + cy * 1.0, -1.0 + cz * 1.0);
                    }
                }
            }

            return gaussianPoints;
        }
    }
}
