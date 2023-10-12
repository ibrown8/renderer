pub mod rasterize;
use geometric_algebra::epga3d::{Point, Plane, Line};//, Line, RegressiveProduct};
use geometric_algebra::{OuterProduct, Transformation, RegressiveProduct};
use core::mem::MaybeUninit;
use sdl2::pixels::Color;
use core::ops::{Fn, Add, Mul, Sub};
use core::mem::take;
//rayon to emulate gpu parallelism
use rayon::iter::{ParallelIterator, IndexedParallelIterator};
type VertexBuffer<V> = Vec<V>;
pub type IndexBuffer = Vec<(u32, u32, u32)>;

#[inline]
pub fn lerp<V>(v0 : V, v1 : V, t : f32) -> V
    where 
        V: Add<Output = V> + Mul<f32, Output = V>
{
      v0 * (1.0 - t) + v1 * t
}

pub trait VertexTrait : Clone + Copy + Send + Add<Output = Self> + Mul<f32, Output = Self> {}
impl VertexTrait for f32 {}
#[repr(C)]
#[derive(Clone, Copy)]
pub struct VertexIn {
    pub vertex_id : u32,
    pub instance_id : u32,
}
#[repr(C)]
#[derive(Clone, Copy)]
pub struct VertexOut {
    pub pos : Point,
    //pub point_size : f32,
    //pub clip_distance : 
}
#[derive(Clone, Copy)]
pub struct PrimPoint<V: VertexTrait>(pub VertexOut, pub V);
//pub type VertexBufferOut<V> = Vec<MaybeUninit<(VertexOut, V)>>;
pub type VertexBufferOut<V: VertexTrait> = Vec<PrimPoint<V>>;
pub type VertexBufferOutRef<'a, V: VertexTrait> = &'a [PrimPoint<V>];
#[derive(Clone, Copy)]
pub struct Triangle<V: VertexTrait> {
    p : [PrimPoint<V>; 3],
}
impl<V: VertexTrait> Triangle<V> {
    pub fn from_points(p0 : PrimPoint<V>, p1 : PrimPoint<V>, p2 : PrimPoint<V>) -> Self {
        Self {p: [p0, p1, p2]}    
    }
}
#[derive(Copy, Clone)]
pub struct ClipPlanes {
    //p : planes in this order top, left, bottom, right
    view_frustum : [Plane; 4],
    near : f32,
    far : f32
}
impl ClipPlanes {
    pub fn from_view_port(width : f32, height : f32, near : f32, far : f32) -> Self {
        let eye = Point::new(1.0, 0.0, 0.0, 0.0);
        let upper_left = Point::new(1.0, width/2.0, height/2.0, far);
        let upper_right = Point::new(1.0, -width/2.0, height/2.0, far);
        let lower_left = Point::new(1.0, width/2.0, -height/2.0, far);
        let lower_right = Point::new(1.0, -width/2.0, -height/2.0, far);
        let top_line = upper_right.regressive_product(upper_left);
        let left_line = upper_left.regressive_product(lower_left);
        let bottom_line = lower_left.regressive_product(lower_right);
        let right_line = lower_right.regressive_product(upper_right);
        let top = top_line.regressive_product(eye);
        let left = left_line.regressive_product(eye);
        let bottom = bottom_line.regressive_product(eye);
        let right = right_line.regressive_product(eye);
        Self {
            view_frustum : [top, left, bottom, right],
            near, far
        } 
    }
}
//https://en.wikipedia.org/wiki/Sutherland%E2%80%93Hodgman_algorithm
impl ClipPlanes {
    pub fn clip_triangle<V: VertexTrait>(&self, tri : &Triangle<V>) -> Vec<PrimPoint<V>> {
        let mut outlist = tri.p.to_vec();
        for (j, plane) in self.view_frustum.iter().enumerate() {
            let inlist = take(&mut outlist);
            for (i, vertex) in inlist.iter().enumerate() {
                let len = inlist.len();
                let prev_vertex = inlist[(i + len - 1) % inlist.len()];
                let line = vertex.0.pos.regressive_product(prev_vertex.0.pos);
                let intersecting_point = line.outer_product(*plane);
                let d0 = plane.regressive_product(vertex.0.pos);
                //TODO: Interpolate other aspects of the vertex
                println!("Distance d0 is {} in iter {},{}", d0, j, i);
                let d1 = plane.regressive_product(prev_vertex.0.pos);
                if d0 <= 0.0 {
                    println!("Distance d1 is {} in iter {},{}", d1, j, i);
                    if d1 > 0.0 {
                        let new_vertex = PrimPoint(VertexOut{pos: intersecting_point}, lerp(vertex.1, prev_vertex.1, (d0)/(d0 + d1)));
                        outlist.push(new_vertex);
                    } else {
                        outlist.push(*vertex);
                    }
                } else if d1 <= 0.0 {
                    outlist.push(*vertex);
                }
            }
        }
        outlist
    }
}
//Captured Variable correspond to uniform
//pub trait VertexShader<V> = Fn(V) + Send + Sync;
//pub trait FragmentShader<F> = Fn(F) -> Color + Send + Sync;
#[inline]
pub fn run_vertex_shader<V0, V1, I, S>(vertex_iter : I, vertices_out : &mut VertexBufferOut<V1>, vertex_shader : S) where
    V0 : Copy + Clone + Send,
    V1 : VertexTrait,
    I : IndexedParallelIterator<Item = V0>,
    S : Fn(V0) -> PrimPoint<V1> + Send + Sync
{
    vertex_iter.map(vertex_shader).collect_into_vec(vertices_out);
}
//Foward facing is counter clockwise
pub fn assemble_triangles<V: VertexTrait>(vertices : &[PrimPoint<V>], indices : &[(u32, u32, u32)]) -> Vec<Triangle<V>> {
    indices.iter().map(|index|{
        let p0 = vertices[index.0 as usize];
        let p1 = vertices[index.1 as usize];
        let p2 = vertices[index.2 as usize];
        Triangle::from_points(p0, p1, p2)
    }).collect()
}
//https://en.wikipedia.org/wiki/Sutherland%E2%80%93Hodgman_algorithm
pub fn perform_clipping<V: VertexTrait>(triangles : &[Triangle<V>], clip_planes : &ClipPlanes) -> Vec<Vec<PrimPoint<V>>> {
    triangles.iter().map(|tri| clip_planes.clip_triangle(tri)).collect()
}
//fn flatten<V: Clone + Copy>
#[inline]
pub fn run_fragment_shader<F, I, S>(frame_buffer : &mut [u8], z_buffer : &mut [f32], fragment_iter : I, fragment_shader : S) where
    I : ParallelIterator<Item = F>,
    S : Fn(F) + Send + Sync
{
    fragment_iter.for_each(fragment_shader);
    /*  frame_buffer.chunks_exact_mut(pitch).enumerate.zip(z_buffer.chuncks_exact_mut()).for_each(|y, scanline, z_line|{
        scanline.chunks_exact_mut(3).enumerate().zip(z_line.parallel_iter_mut()).for_each(|x, pixel, z_index|{
            for triangle in Triangles {
                //check if pixel is in trianngle
            }      
        }
    }); */
    
}
pub struct FragmentIn {
    pub frag_coordinate : [f32; 4],
    pub front_facing : bool,
    pub point_coordinate : [f32; 2]
}

/* pub fn transform_points(motor : &Motor, points_in : &[Vertex], points_out : &mut [Vertex]){
    for (p_in, p_out) in points_in.iter().zip(points_out.iter_mut()){
        *p_out = motor.transformation(*p_in);
    }
} */

