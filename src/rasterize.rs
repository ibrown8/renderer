use crate::{Triangle, PrimPoint, VertexTrait};
use core::cmp::{min, max};
use sdl2::pixels::Color;
struct Point2 {
    x : f32,
    y : f32
}
struct BoundingBox {
    upper : Point2, //Upper Left
    lower : Point2 //Bottom Right
}
fn clamp(val : f32, lower : f32, upper : f32) -> f32 {
    lower.max(val.min(upper))
}
//https://erkaman.github.io/posts/fast_triangle_rasterization.html
#[inline]
fn line_test(v0 : &Point2, v1 : &Point2, p : &Point2) -> bool {
    (((v1.x - v0.x) * (p.y - v0.y)) - ((v1.y - v0.y) * (p.x - v0.x))) >= 0.0
}
//i16 instead of usize to make conversion to and float easier
pub fn rasterize_triangle<V: VertexTrait>(tri : &Triangle<V>, frame_width : i16, frame_height : i16, frame_buffer : &mut [u8], pitch : usize, color : Color){
    let v0 = tri.p[0].0.pos.group0();
    let v1 = tri.p[1].0.pos.group0();
    let v2 = tri.p[2].0.pos.group0();
    let v0 = Point2{x : -v0[1], y : v0[2]};
    let v1 = Point2{x : -v1[1], y : v1[2]};
    let v2 = Point2{x : -v2[1], y : v2[2]};
    let x_min = clamp(v0.x.min(v1.x.min(v2.x)), -1.0, 1.0);
    let y_min = clamp(v0.y.min(v1.y.min(v2.y)), -1.0, 1.0);
    let x_max = clamp(v0.x.max(v1.x.max(v2.x)), -1.0, 1.0);
    let y_max = clamp(v0.y.max(v1.y.max(v2.x)), -1.0, 1.0);
    let bounding_box = BoundingBox {
        upper : Point2{x : x_min, y : y_max},
        lower : Point2{x : x_max, y : y_min}
    };
    let screen_center = Point2{x : frame_width as f32/2.0, y : frame_height as f32/2.0};
    let y0 : usize = unsafe {(screen_center.y - (bounding_box.upper.y * screen_center.y)).to_int_unchecked() };//should be positive!
    let y1 : usize = unsafe {(screen_center.y - (bounding_box.lower.y * screen_center.y)).to_int_unchecked() };//should be positive!
    let x0 : usize = unsafe {(screen_center.x + (bounding_box.upper.x * screen_center.x)).to_int_unchecked() };//should be positive!
    let x1 : usize = unsafe {(screen_center.x + (bounding_box.upper.x * screen_center.x)).to_int_unchecked() };//should be positive!
    println!("rasterizing in range y : {}-{}, x : {}-{}", y0, y1, x0, x1);
    for y in y0..y1 {
        let row = &mut frame_buffer[(y * pitch)..];
        for x in y0..y1 {
            let offset = x * 3;
            row[offset + 0] = color.r;
            row[offset + 1] = color.g;
            row[offset + 2] = color.b;
        }
    }
}
