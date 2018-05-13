﻿namespace WpfApp1
{
    using System.Collections.Generic;
    using System.Windows.Media;
    using System.Windows.Media.Media3D;

    using HelixToolkit.Wpf;

    /// <summary>
    /// Provides a ViewModel for the Main window.
    /// </summary>
    public class MainViewModel
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="MainViewModel"/> class.
        /// </summary>
        public MainViewModel()
        {
            // Create a model group
            //var modelGroup = new Model3DGroup();

            //// Create a mesh builder and add a box to it
            //var meshBuilder = new MeshBuilder();
            ////meshBuilder.AddBox(new Point3D(0, 0, 1), 1, 2, 0.5);
            ////meshBuilder.AddBox(new Rect3D(0, 0, 1.2, 0.5, 1, 0.4));
            ////meshBuilder.AddEdges(new List<Point3D>() { new Point3D(0, 0, 0), new Point3D(3, 3, 3), new Point3D(6, 6, 6) }, new List<int>() { 0, 1, 2 }, 2, 1);
            //meshBuilder.AddPipe(new Point3D(0, 0, 0), new Point3D(3, 3, 3), 1, 10, 0);
            ////meshBuilder.add

            //// Create a mesh from the builder (and freeze it)
            //var mesh = meshBuilder.ToMesh();

            //// Create some materials
            //var greenMaterial = MaterialHelper.CreateMaterial(Colors.Green);
            ////var redMaterial = MaterialHelper.CreateMaterial(Colors.Red);
            ////var blueMaterial = MaterialHelper.CreateMaterial(Colors.Blue);
            //var insideMaterial = MaterialHelper.CreateMaterial(Colors.Yellow);

            //// Add 3 models to the group (using the same mesh, that's why we had to freeze it)
            ////modelGroup.Children.Add(new GeometryModel3D { Geometry = mesh, Material = greenMaterial, BackMaterial = insideMaterial });
            ////modelGroup.Children.Add(new GeometryModel3D { Geometry = mesh, Transform = new TranslateTransform3D(-2, 0, 0), Material = redMaterial, BackMaterial = insideMaterial });
            ////modelGroup.Children.Add(new GeometryModel3D { Geometry = mesh, Transform = new TranslateTransform3D(2, 0, 0), Material = blueMaterial, BackMaterial = insideMaterial });
            //modelGroup.Children.Add(new GeometryModel3D { Geometry = mesh, Material = greenMaterial, BackMaterial = insideMaterial });

            //// Set the property, which will be bound to the Content property of the ModelVisual3D (see MainWindow.xaml)
            //this.Model = modelGroup;

            
        }

        /// <summary>
        /// Gets or sets the model.
        /// </summary>
        /// <value>The model.</value>
        public Model3D Model { get; set; }
    }
}