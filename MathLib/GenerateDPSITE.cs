using System;

namespace MathLib
{
    public class GenerateDPSITE
    {
        public static Tuple<int, int>[] niti = GetNiTi();

        public static Func<double, double, double>[] PHIs = new Func<double, double, double>[]
        {
            (n, t) => 0.25 * (1 + n * niti[0].Item1) * (1 + t * niti[0].Item2) * (n * niti[0].Item1 + t * niti[0].Item2 - 1),
            (n, t) => 0.25 * (1 + n * niti[1].Item1) * (1 + t * niti[1].Item2) * (n * niti[1].Item1 + t * niti[1].Item2 - 1),
            (n, t) => 0.25 * (1 + n * niti[2].Item1) * (1 + t * niti[2].Item2) * (n * niti[2].Item1 + t * niti[2].Item2 - 1),
            (n, t) => 0.25 * (1 + n * niti[3].Item1) * (1 + t * niti[3].Item2) * (n * niti[3].Item1 + t * niti[3].Item2 - 1),
            (n, t) => 0.5 * (1 - n * n) * (1 + t * niti[4].Item2),
            (n, t) => 0.5 * (1 - t * t) * (1 + n * niti[5].Item1),
            (n, t) => 0.5 * (1 - n * n) * (1 + t * niti[6].Item2),
            (n, t) => 0.5 * (1 - t * t) * (1 + n * niti[7].Item1),
        };

        public static double[,,] Generate()
        {
            var gauss = GetGaussNode();
            var res = new double[9, 2, 8];
            for (int i = 0; i < 9; i++)
            {
                for (int j = 0; j < 8; j++)
                {
                    res[i, 0, j] = DN(j, gauss[i, 0], gauss[i, 1]);
                    res[i, 1, j] = DT(j, gauss[i, 0], gauss[i, 1]);
                }
            }
            return res;
        }

        public static double DN(int j, double p1, double p2)
        {
            var nt = GetNiTi();
            if (j < 4)
            {
                return 0.25 * nt[j].Item1 * (nt[j].Item2 * p2 + 1) * (2 * nt[j].Item1 * p1 + nt[j].Item2 * p2);
            }
            else if (j == 4 || j == 6)
            {
                return -p1 * (nt[j].Item2 * p2 + 1);
            }
            else
            {
                return 0.5 * nt[j].Item1 * (1 - p2 * p2);
            }
        }

        public static double DT(int j, double p1, double p2)
        {
            var nt = GetNiTi();
            var gauss = GetGaussNode();
            if (j < 4)
            {
                return 0.25 * nt[j].Item2 * (nt[j].Item1 * p1 + 1) * (2 * nt[j].Item2 * p2 + nt[j].Item1 * p1);
            }
            else if (j == 4 || j == 6)
            {
                return 0.5 * nt[j].Item2 * (1 - p1 * p1);
            }
            else
            {
                return -p2 * (nt[j].Item1 * p1 + 1);
            }
        }

        // TODO: return fking multi dimensional arrray
        public static Tuple<int, int>[] GetNiTi()
        {
            var result = new Tuple<int, int>[8];
            result[0] = new Tuple<int, int>(-1, -1);
            result[1] = new Tuple<int, int>(1, -1);
            result[2] = new Tuple<int, int>(1, -1);
            result[3] = new Tuple<int, int>(-1, 1);
            result[4] = new Tuple<int, int>(0, -1);
            result[5] = new Tuple<int, int>(1, 0);
            result[6] = new Tuple<int, int>(0, 1);
            result[7] = new Tuple<int, int>(-1, 0);
            return result;
        }

        public static double[,] GetGaussNode()
        {
            double[,] result = new double[9,2];
            var sq = Math.Sqrt(0.6);
            double[] x123 = new double[] { -sq, 0, sq };
            int ij = 0;
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    result[ij, 0] = x123[i];
                    result[ij, 1] = x123[j];
                    ++ij;
                }
            }
            return result;
        }
    }
}
