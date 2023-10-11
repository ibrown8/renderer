use geometric_algebra::epga3d::{Point, Plane, Line};//, Line, RegressiveProduct};
use geometric_algebra::{OuterProduct, Transformation, RegressiveProduct};
use core::mem::MaybeUninit;
use sdl2::pixels::Color;
use core::ops::Fn;
use core::mem::take;
//rayon to emulate gpu parallelism
use rayon::iter::{ParallelIterator, IndexedParallelIterator};
type VertexBuffer<V> = Vec<V>;
type IndexBuffer = Vec<(u32, u32, u32)>;

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
pub struct PrimPoint<V: Clone + Copy + Send>(VertexOut, V);
//pub type VertexBufferOut<V> = Vec<MaybeUninit<(VertexOut, V)>>;
pub type VertexBufferOut2<V> = Vec<PrimPoint<V>>;
pub type VertexBufferOut2Ref<'a, V> = &'a [PrimPoint<V>];
#[derive(Clone, Copy)]
pub struct Triangle<V: Clone + Copy + Send> {
    p : [PrimPoint<V>; 3],
}
impl<V: Clone + Copy + Send> Triangle<V> {
    pub fn from_points(p0 : PrimPoint<V>, p1 : PrimPoint<V>, p2 : PrimPoint<V>) -> Self {
        Self {p: [p0, p1, p2]}    
    }
}
pub struct ClipPlanes {
    //p : planes in this order top, left, bottom, right
    pub view_frustum : [Plane; 4],
    pub near : f32,
    pub far : f32
}
//https://en.wikipedia.org/wiki/Sutherland%E2%80%93Hodgman_algorithm
impl ClipPlanes {
    pub fn clip_triangle<V: Clone + Copy + Send>(&self, tri : &Triangle<V>) -> Vec<PrimPoint<V>> {
        let mut outlist = tri.p.to_vec();
        for (j, plane) in self.view_frustum.iter().enumerate() {
            let inlist = take(&mut outlist);
            for (i, vertex) in inlist.iter().enumerate() {
                let prev_vertex = inlist[(i - 1) % inlist.len()];
                let line = vertex.0.pos.regressive_product(prev_vertex.0.pos);
                let intersecting_point = line.outer_product(*plane);
                let d0 = plane.regressive_product(vertex.0.pos);
                //TODO: Interpolate other aspects of the vertex
                println!("Distance d0 is {} in iter {}{}", d0, j, i);
                let d1 = plane.regressive_product(prev_vertex.0.pos);
                if d0 <= 0.0 {
                    println!("Distance d1 is {} in iter {}{}", d1, j, i);
                    if d1 > 0.0 {
                        outlist.push(intersecting_point);
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

/* #[inline]
pub fn run_vertex_shader<V, I, S>(vertex_iter : I, vertex_shader : S) where
    I : ParallelIterator<Item = V>,
    S : Fn(V) + Send + Sync
{
    vertex_iter.for_each(vertex_shader);
} */

#[inline]
pub fn run_vertex_shader2<V0, V1, I, S>(vertex_iter : I, vertex_shader : S, vertices_out : &mut VertexBufferOut2<V1>) where
    V0 : Copy + Clone + Send,
    V1 : Copy + Clone + Send,
    I : IndexedParallelIterator<Item = V0>,
    S : Fn(V0) -> PrimPoint<V1> + Send + Sync
{
    vertex_iter.map(vertex_shader).collect_into_vec(vertices_out);
}
//Foward facing is counter clockwise
pub fn assemble_triangles<V: Clone + Copy + Send>(vertices : &[PrimPoint<V>], indices : &[(u32, u32, u32)]) -> Vec<Triangle<V>> {
    indices.iter().map(|index|{
        let p0 = vertices[index.0 as usize];
        let p1 = vertices[index.1 as usize];
        let p2 = vertices[index.2 as usize];
        Triangle::from_points(p0, p1, p2)
    }).collect()
}
//https://en.wikipedia.org/wiki/Sutherland%E2%80%93Hodgman_algorithm
pub fn perform_clipping<V: Clone + Copy + Send>(triangles : &[Triangle<V>], clip_planes : &ClipPlanes) -> Vec<Vec<PrimPoint<V>>> {
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
//Bresenhams algorithm
//https://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm
fn draw_line_low<P>(x0 : i32, y0 : i32, x1 : i32, y1 : i32, plot : P)
    where P : FnMut(i32, i32)
{
    let dx = x1 - x0;
    let (dy, yi) = if y1 < y0 {
        (y0 - y1, 1)
    } else {
        (y1 - y0, -1)
    };
    let mut d = (2 * dy) - dx;
    let mut y = y0;
    for x in x0..x1 {
        plot(x, y); 
        if d > 0 {
            y += yi;
            d += 2 * (dy - dx);
        } else {
            d += 2 * dy;
        }
    }
}
fn draw_line_high<P>(x0 : i32, y0 : i32, x1 : i32, y1 : i32, plot : P)
    where P : FnMut(i32, i32)
{
    let dy = y1 - y0;
    let (dx, xi) = if x1 < x0 {
        (x0 - x1, 1)
    } else {
        (x1 - x0, -1)
    };
    let mut d = (2 * dx) - dy;
    let mut x = x0;
    for y in y0..y1 {
        plot(x, y); 
        if d > 0 {
            x += xi;
            d += 2 * (dx - dy);
        } else {
            d += 2 * dx;
        }
    }
}
pub fn draw_line<P>(x0 : i32, y0 : i32, x1 : i32, y1 : i32, plot : P)
    where P : FnMut(i32, i32)
{
    if (y1 - y0).abs() > (x1 - x0).abs() {
        if x0 > x1 {
            draw_line_low(x1, y1, x0, y0, plot);
        } else {
            draw_line_low(x0, y0, x1, y1, plot);
        }
    } else {
        if y0 > y1 {
            draw_line_high(x1, y1, x0, y0, plot);
        } else {
            draw_line_high(x0, y0, x1, y1, plot);
        }
    }
}

/* pub fn rasterize_traingles(points : &[Point], indices : &[(u32, u32, u32)], buffer : &mut [u8], pitch : u32, height : u32, width : u32, color : Color){
    let center_x = width / 2;
    let center_y = height / 2;
    for index in indices {
        let p0 = <[f32; 4]>::from(points[index.0 as usize]);
        let p1 = <[f32; 4]>::from(points[index.1 as usize]);
        let p2 = <[f32; 4]>::from(points[index.2 as usize]);
        let p0x = p0[0] as i32;
        let p0y = p0[1] as i32;
        let p1x = p1[0] as i32;
        let p1y = p1[1] as i32;
        let p2x = p2[0] as i32;
        let p2y = p2[1] as i32;
        let plot = |x : i32, y : i32| {
            //positive x has pixel indices
            let screen_x : usize = center_x - x;
            let screen_y : usize = center_y - y;
            let offset = screen_y * pitch + screen_x * 3;
            buffer[offset + 0] = color.r;
            buffer[offset + 1] = color.g;
            buffer[offset + 2] = color.b;
        };
        draw_line(p0x, p0y, p1x, p1y, plot);
        draw_line(p1x, p1y, p2x, p2y, plot);
        draw_line(p2x, p2y, p0x, p0y, plot);
    }
} */
    
