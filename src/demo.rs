use renderer::{run_vertex_shader, run_fragment_shader, VertexOut};
use geometric_algebra::{Point}
use core::mem::MaybeUninit;
pub fn main(){
    let sdl_context = sdl2::init().unwrap();
    let video_subsystem = sdl_context.video().unwrap();
    let window = video_subsystem.window("renderer_demo", 640, 480).
        position_centered().
        build().
        unwrap();
    let mut canvas = window.into_canvas().build().unwrap();
    println!("Canvas Created");
    let texture_builder = canvas.texture_creator();
    let points = 
    let mut texture = texture_builder
        .create_texture_streaming(PixelFormatEnum::RGB24, 640, 480)
        .map_err(|e| e.to_string())?;
    let vertices = vec![
        Point::new(1.0, 0.5, -0.5, 0.0),
        Point::new(1.0, -0.5, -0.5, 0.0),
        Point::new(1.0, 0.0, 0.5, 0.0)
    ];
    let mut vertices_out : VertexBufferOut<()>
    let indices : IndexBuffer = vec![(0, 1, 2)];
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
            let vertices_out = run_vertex_shader2(vertex_iter, |vertex_in|{
                PrimPoint(VertexOut{p : vertex_in}, ())
            })
            let triangles = assemble_triangles(&vertices_out, &indices);
            let triangles2 = perform_clipping(&triangles);
        }
    }
}
