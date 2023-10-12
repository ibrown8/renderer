use renderer::{run_vertex_shader, run_fragment_shader, VertexOut, VertexTrait};
use renderer::{PrimPoint, VertexBufferOut, IndexBuffer, assemble_triangles, perform_clipping};
use renderer::rasterize::rasterize_triangle;
use geometric_algebra::epga3d::{Point, Plane};
use renderer::ClipPlanes;
use core::mem::MaybeUninit;
use sdl2::event::Event;
use sdl2::keyboard::Keycode;
use sdl2::pixels::{PixelFormatEnum, Color};
use rayon::iter::{ParallelIterator, IntoParallelRefIterator};
use std::time::Duration;
pub fn main(){
    //Eye is located at 0.0, 0.0, 0.0, 0.0
    let sdl_context = sdl2::init().unwrap();
    let video_subsystem = sdl_context.video().unwrap();
    let window = video_subsystem.window("renderer_demo", 640, 480).
        position_centered().
        build().
        unwrap();
    let mut canvas = window.into_canvas().build().unwrap();
    canvas.set_draw_color(Color::RGB(0, 255, 255));
    let mut event_pump = sdl_context.event_pump().unwrap();
    println!("Canvas Created");
    let texture_builder = canvas.texture_creator();
    let mut texture = texture_builder.create_texture_streaming(PixelFormatEnum::RGB24, 640, 480).unwrap();
    let vertices = vec![
        Point::new(1.0, 0.5, -0.5, -1.0),
        Point::new(1.0, -0.5, -0.5, -1.0),
        Point::new(1.0, 0.0, 0.5, -1.0)
    ];
    //let clip_planes = ClipPlanes::from_view_port(1.3333, 1.000, 0.5, 6.0);
    let mut vertices_out : VertexBufferOut<f32> = VertexBufferOut::new();
    let indices : IndexBuffer = vec![(0, 1, 2)];
    'running : loop {
        for event in event_pump.poll_iter(){
            match event {
                Event::Quit{..} | 
                Event::KeyDown {keycode : Some(Keycode::Escape), ..} => {
                    println!("Loop Terminating");
                    break 'running;
                }
                _ => {}
            }
            texture.with_lock(None, |buffer: &mut [u8], pitch: usize| {
                //let vertex_iter = vertices.par_iter().zip(vertices_out.par_iter_mut());
                for y in 0..480 {
                    let row = &mut buffer[(y * pitch)..];
                    for x in 0..1920 {
                        buffer[x] = 0;
                    }
                }
                run_vertex_shader(vertices.par_iter(), &mut vertices_out, |vertex_in|{
                    PrimPoint(VertexOut{pos : *vertex_in}, 0.0)
                });
                let triangles = assemble_triangles(&vertices_out, &indices);
                //let triangles2 = perform_clipping(&triangles, &clip_planes);
                for tri in &triangles {
                    rasterize_triangle(tri, 640, 480, buffer, pitch, Color::RGB(255, 0, 0));
                }
            });
            canvas.clear();
            canvas.copy(&texture, None, None).unwrap();
            canvas.present();
            ::std::thread::sleep(Duration::new(0, 1_000_000_000u32 / 30));
        }
    }
}
