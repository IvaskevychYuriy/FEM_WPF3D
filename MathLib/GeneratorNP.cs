using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MathLib
{
    public class GeneratorNP
    {
        public static double[,] Generate(int[,,] globalToAKT,  int nx, int ny, int nz)
        {
            var nt = new double[nx * ny * nz, 20];
            int counter = 0;
            for (int iz = 0; iz < nz; iz++)
            {
                for (int iy = 0; iy < ny; iy++)
                {
                    for (int ix = 0; ix < nx; ix++)
                    {
                        nt[counter, 0] = globalToAKT[2 * ix, 2 * iy, 2 * iz];
                        nt[counter, 1] = globalToAKT[2 * ix + 1, 2 * iy, 2 * iz];
                        nt[counter, 2] = globalToAKT[2 * ix + 2, 2 * iy, 2 * iz];
                        nt[counter, 3] = globalToAKT[2 * ix, 2 * iy + 1, 2 * iz];
                        nt[counter, 4] = globalToAKT[2 * ix + 2, 2 * iy + 1, 2 * iz];
                        nt[counter, 5] = globalToAKT[2 * ix, 2 * iy + 2, 2 * iz];
                        nt[counter, 6] = globalToAKT[2 * ix + 1, 2 * iy + 2, 2 * iz];
                        nt[counter, 7] = globalToAKT[2 * ix + 2, 2 * iy + 2, 2 * iz];

                        nt[counter, 8] = globalToAKT[2 * ix, 2 * iy, 2 * iz + 1];
                        nt[counter, 9] = globalToAKT[2 * ix + 2, 2 * iy, 2 * iz + 1];
                        nt[counter, 10] = globalToAKT[2 * ix, 2 * iy + 2, 2 * iz + 1];
                        nt[counter, 11] = globalToAKT[2 * ix + 2, 2 * iy + 2, 2 * iz + 1];

                        nt[counter, 12] = globalToAKT[2 * ix, 2 * iy, 2 * iz + 2];
                        nt[counter, 13] = globalToAKT[2 * ix + 1, 2 * iy, 2 * iz + 2];
                        nt[counter, 14] = globalToAKT[2 * ix + 2, 2 * iy, 2 * iz + 2];
                        nt[counter, 15] = globalToAKT[2 * ix, 2 * iy + 1, 2 * iz + 2];
                        nt[counter, 16] = globalToAKT[2 * ix + 2, 2 * iy + 1, 2 * iz + 2];
                        nt[counter, 17] = globalToAKT[2 * ix, 2 * iy + 2, 2 * iz + 2];
                        nt[counter, 18] = globalToAKT[2 * ix + 1, 2 * iy + 2, 2 * iz + 2];
                        nt[counter, 19] = globalToAKT[2 * ix + 2, 2 * iy + 2, 2 * iz + 2];
                        counter++;
                    }
                }
            }
            return nt;
        }
    }
}
