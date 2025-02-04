
import dearpygui.dearpygui as dpg
from math import sin

dpg.create_context()

# creating data
datax = []
datay = []
for i in range(0, 500):
    datax.append(i / 1000)
    datay.append(0.5 + 0.5 * sin(50 * i / 1000))

with dpg.window(label="Tutorial"):
    # create plot
    with dpg.plot(label="Line Series", height=400, width=400):
        # optionally create legend
        dpg.add_plot_legend()

        # REQUIRED: create x and y axes
        dpg.add_plot_axis(dpg.mvXAxis, label="x")
        dpg.add_plot_axis(dpg.mvYAxis, label="y", tag="y_axis")

        # series belong to a y axis
        dpg.add_line_series(datax, datay, label="0.5 + 0.5 * sin(x)", parent="y_axis")

dpg.create_viewport(title='Custom Title', width=800, height=600)
dpg.setup_dearpygui()
dpg.show_viewport()
dpg.start_dearpygui()
dpg.destroy_context()

