using System;
using System.Collections.Generic;

namespace MathLib
{
    public class GenerateDPSITE
    {
        public static Dictionary<int, int[]> gauss = new Dictionary<int, int[]> {
            {0,  new int[]{0,0,0} },
            {1,  new int[]{2,0,0} },
            {2,  new int[]{2,2,0} },
            {3,  new int[]{0,2,0} },

            {4,  new int[]{0,0,2} },
            {5,  new int[]{2,0,2} },
            {6,  new int[]{2,2,2} },
            {7,  new int[]{0,2,2} },

            {8,  new int[]{1,0,0} },
            {9,  new int[]{2,1,0} },
            {10, new int[]{1,2,0} },
            {11, new int[]{0,1,0} },

            {12, new int[]{0,0,1} },
            {13, new int[]{2,0,1} },
            {14, new int[]{2,2,1} },
            {15, new int[]{0,2,1} },

            {16, new int[]{1,0,2} },
            {17, new int[]{2,1,2} },
            {18, new int[]{1,2,2} },
            {19, new int[]{0,1,2} },
        };

        public static double[,,] Generate()
        {
            var res = new double[9, 2, 8];
            for (int i = 0; i < 9; i++)
            {
                for (int j = 0; j < 8; j++)
                {
                    res[i, 0, j] = DN(i, j);
                    res[i, 1, j] = DT(i, j);
                }
            }
            return res;
        }

        public static double DN(int i, int j)
        {
            var nt = GetNiTi();
            if (j < 4)
            {
                return 1 / 4 * nt[j].Item1 * (nt[j].Item2 * gauss[i][1] + 1) * (2 * nt[j].Item1 * gauss[i][0] + nt[j].Item2 * gauss[i][1]);
            }
            else if (j == 5 || j == 7)
            {
                return nt[j].Item2 * gauss[i][1] * gauss[i][0] + gauss[i][0];
            }
            else
            {
                return 1 / 2 * nt[j].Item1 * (gauss[i][1] * gauss[i][1] + 1);
            }
        }

        public static double DT(int i, int j)
        {
            var nt = GetNiTi();
            if (j < 4)
            {
                return 1 / 4 * nt[j].Item2 * (nt[j].Item1 * gauss[i][0] + 1) * (2 * nt[j].Item2 * gauss[i][1] + nt[j].Item1 * gauss[i][0]);
            }
            else if (j == 5 || j == 7)
            {
                return 1 / 2 * nt[j].Item2 * (gauss[i][0] * gauss[i][0] + 1);
            }
            else
            {
                return nt[j].Item1 * gauss[i][1] * gauss[i][0] + gauss[i][1];
            }
        }

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
