using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using MathLib;
using Microsoft.AspNetCore.Mvc;

namespace Backend.Controllers
{
    [Route("api/result")]
    public class ValuesController : Controller
    {
        // GET api/values
        [HttpGet]
        public object Get()
        {
            MainSolution solve = new MainSolution();
            var result = solve.Start();
            return new
            {
                start = result.Item1,
                end = result.Item2
            };
        }
    }
}
