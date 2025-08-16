from PIL import Image, ImageDraw, ImageFont
import os


def create_checkerboard(rows, cols, square_size, colors, texts, output_file):
    width = cols * square_size
    height = rows * square_size
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)

    font = ImageFont.truetype("Arial.ttf", 35)

    for row in range(rows):
        for col in range(cols):
            color = colors[(row * cols + (col + row)) % len(colors)]

            x0 = col * square_size
            y0 = row * square_size
            x1 = x0 + square_size
            y1 = y0 + square_size

            draw.rectangle([x0, y0, x1, y1], fill=color)

            text = texts[row][col]
            text_x = x0 + (square_size) // 3
            text_y = y0 + (square_size) // 3
            draw.text((text_x, text_y), text, fill="black", font=font)

    img.save(os.path.join("resources", output_file))
    print(f"Checkerboard saved as {output_file}")


rows = 8
cols = 8
square_size = 1024 // 8
texts = [[f"{chr(65 + row)}{col + 1}" for col in range(cols)] for row in range(rows)]

output_file = "celestial.png"
colors = ["#E8A0B6", "#96B5D6", "#F2E394", "#88C099"]
create_checkerboard(rows, cols, square_size, colors, texts, output_file)

output_file = "sphere.png"
colors = ["#FFB5A7", "#FCD5CE", "#F8E5A8", "#B7CFA2"]
create_checkerboard(rows, cols, square_size, colors, texts, output_file)

output_file = "disk.png"
colors = ["#FF6B35", "#F7931E", "#C5282F", "#FF9F1C"]
create_checkerboard(rows, cols, square_size, colors, texts, output_file)
