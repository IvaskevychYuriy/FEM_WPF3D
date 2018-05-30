using System;
using System.Collections.Generic;
using System.Text;

namespace MathLibrary
{
    public class Point
    {
        Point(int x, int y, int z)
        {
            X = x;
            Y = y;
            Z = z;
        }

        public int X { get; set; }
        public int Y { get; set; }

        public int Z { get; set; }
    }
}
