namespace MathLib
{
    public class GeneratorNP
    {
        public static int[,] Generate(int[,,] globalToAKT,  int nx, int ny, int nz)
        {
            var nt = new int[20, nx * ny * nz];
            int counter = 0;
            for (int iz = 0; iz < nz; iz++)
            {
                for (int iy = 0; iy < ny; iy++)
                {
                    for (int ix = 0; ix < nx; ix++)
                    {
                        nt[0, counter] = globalToAKT[2 * ix, 2 * iy, 2 * iz];
                        nt[1, counter] = globalToAKT[2 * ix + 1, 2 * iy, 2 * iz];
                        nt[2, counter] = globalToAKT[2 * ix + 2, 2 * iy, 2 * iz];
                        nt[3, counter] = globalToAKT[2 * ix, 2 * iy + 1, 2 * iz];
                        nt[4, counter] = globalToAKT[2 * ix + 2, 2 * iy + 1, 2 * iz];
                        nt[5, counter] = globalToAKT[2 * ix, 2 * iy + 2, 2 * iz];
                        nt[6, counter] = globalToAKT[2 * ix + 1, 2 * iy + 2, 2 * iz];
                        nt[7, counter] = globalToAKT[2 * ix + 2, 2 * iy + 2, 2 * iz];
                             
                        nt[8, counter] = globalToAKT[2 * ix, 2 * iy, 2 * iz + 1];
                        nt[9, counter] = globalToAKT[2 * ix + 2, 2 * iy, 2 * iz + 1];
                        nt[10, counter] = globalToAKT[2 * ix, 2 * iy + 2, 2 * iz + 1];
                        nt[11, counter] = globalToAKT[2 * ix + 2, 2 * iy + 2, 2 * iz + 1];
                             
                        nt[12, counter] = globalToAKT[2 * ix, 2 * iy, 2 * iz + 2];
                        nt[13, counter] = globalToAKT[2 * ix + 1, 2 * iy, 2 * iz + 2];
                        nt[14, counter] = globalToAKT[2 * ix + 2, 2 * iy, 2 * iz + 2];
                        nt[15, counter] = globalToAKT[2 * ix, 2 * iy + 1, 2 * iz + 2];
                        nt[16, counter] = globalToAKT[2 * ix + 2, 2 * iy + 1, 2 * iz + 2];
                        nt[17, counter] = globalToAKT[2 * ix, 2 * iy + 2, 2 * iz + 2];
                        nt[18, counter] = globalToAKT[2 * ix + 1, 2 * iy + 2, 2 * iz + 2];
                        nt[19, counter] = globalToAKT[2 * ix + 2, 2 * iy + 2, 2 * iz + 2];
                        counter++;
                    }
                }
            }
            return nt;
        }
    }
}
