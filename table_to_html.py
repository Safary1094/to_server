import os
from math import floor
from os.path import join

import pandas as pd

from ngs_utils.utils import mean
from ngs_reporting.rnaseq.style.style_css import table_css_string, page_css_string

# Color heat map
BLUE_HUE = 240
BLUE_OUTER_BRT = 55
BLUE_INNER_BRT = 65

GREEN_HUE = 120
GREEN_OUTER_BRT = 50
GREEN_INNER_BRT = 60

RED_HUE = 0
RED_OUTER_BRT = 50
RED_INNER_BRT = 60

MIN_NORMAL_BRT = 80
MEDIAN_BRT = 100  # just white.


def calc_row_stats(row):
    numbers = sorted([v for v in row],
                     key=lambda a: a if a is not None else -1)  # None is always less than anything
    l = len(numbers)

    row_min = numbers[0]
    row_max = numbers[l - 1]
    row_is_all_values_equal = row_min == row_max
    row_med = numbers[(l - 1) // 2] if l % 2 != 0 else mean([numbers[l // 2], numbers[(l // 2) - 1]])
    q1 = numbers[int(floor((l - 1) // 4))]
    q3 = numbers[int(floor((l - 1) * 3 // 4))]

    d = q3 - q1
    row_low_outer_fence = q1 - 3 * d
    row_low_inner_fence = q1 - 1.5 * d
    row_top_inner_fence = q3 + 1.5 * d
    row_top_outer_fence = q3 + 3   * d

    return row_min, row_max, row_med, row_is_all_values_equal, row_low_outer_fence, row_low_inner_fence, row_top_outer_fence, row_top_inner_fence


def set_cell_colors(row):
    if row.ndim > 1:
        return ''

    #row = data[2:]

    row_min, row_max, row_med, is_all_values_equal, low_outer_fence, low_inner_fence, top_outer_fence, top_inner_fence = calc_row_stats(row)

    # Color heatmap
    bg_color_attr  = []
    #txt_color_attr = []
    for i, r in enumerate(row):
        #if i < 2:
        #    txt_color_attr.append('')
        #    bg_color_attr.append('')
        #    continue

        #text_color = 'black'
        color = 'white'
        if r is not None:
            [top_hue, inner_top_brt, outer_top_brt] = [BLUE_HUE, BLUE_INNER_BRT, BLUE_OUTER_BRT]
            [low_hue, inner_low_brt, outer_low_brt] = [RED_HUE, RED_INNER_BRT, RED_OUTER_BRT]

            if not is_all_values_equal:
                #text_color = 'black'

                # Low outliers
                if r < low_outer_fence and r < row_med:
                    color = get_color(low_hue, outer_low_brt)
                    #text_color = 'white'

                elif r < low_inner_fence and r < row_med:
                    color = get_color(low_hue, inner_low_brt)

                # Normal values
                elif r < row_med:
                    try:
                        k = float(MEDIAN_BRT - MIN_NORMAL_BRT) / (row_med - low_inner_fence)
                    except:
                        pass
                    else:
                        brt = round(MEDIAN_BRT - (row_med - r) * k)
                        color = get_color(low_hue, brt)

                # High outliers
                elif r > top_inner_fence and r > row_med:
                    color = get_color(top_hue, inner_top_brt)

                elif r > top_outer_fence and r > row_med:
                    color = get_color(top_hue, outer_top_brt)
                    #text_color = 'white'

                elif r > row_med:
                    k = float(MEDIAN_BRT - MIN_NORMAL_BRT) / (top_inner_fence - row_med)
                    brt = round(MEDIAN_BRT - (r - row_med) * k)
                    color = get_color(top_hue, brt)

        #txt_color_attr.append('color: {}'.format(text_color))
        bg_color_attr.append('background-color: {}'.format(color))

    return bg_color_attr


def hue2rgb(p, q, t):
    if t < 0: t += 1
    if t > 1: t -= 1
    if t < 1./6: return p + (q - p) * 6 * t
    if t < 1./2: return q
    if t < 2./3: return p + (q - p) * (2./3 - t) * 6
    return p


def hsl2rgb(h, s, l):
    r, g, b = None, None, None

    if s == 0:
        r = g = b = l  # achromatic
    else:
        q = l * (1 + s) if l < 0.5 else l + s - l * s
        p = 2 * l - q
        r = hue2rgb(p, q, h + 1./3)
        g = hue2rgb(p, q, h)
        b = hue2rgb(p, q, h - 1./3)

    return map(int, [round(r * 255), round(g * 255), round(b * 255)])


def get_color(hue, lightness):
    lightness = lightness or 92
    rgb = hsl2rgb(float(hue) / 360, 0.8, float(lightness) / 100)
    hex_rgb = [hex(c)[2:] for c in rgb]
    return '#' + ''.join(hex_rgb)


def table_to_html(table, gradient_cols, title, path):
    table_id = 'Level'
    # gradient_cols = ["ILS38024-PT1-DS1_S1","Karpas299_S4","RT4_S2","RT112_S3"]
    N = len(gradient_cols)
    num_of_cols = len(table.columns)

    str_col_width = "30px"
    max_name = len(max(gradient_cols, key=len))
    str_col_height = str(max_name * 10) + "px"

    styles = [
        dict(selector="td", props=[("padding", "4px")]),
        dict(selector="thead th:first-child", props=[("display", "none")]),
        dict(selector="tbody th:first-child", props=[("display", "none")]),
        dict(selector="thead th", props=[("height", str_col_height)]),
        dict(selector="th:nth-child(-n+" + str(num_of_cols - N + 1) + ")", props=[
                                                    ("position", "relative"),
                                                    ("vertical-align", "bottom"),
                                                    ("text-align", "left")
                                                  ]),
        dict(selector="th:nth-child(n+" + str(num_of_cols - N + 2) + ")", props=[("text-align", "center"),
                                                        ("position", "relative"),
                                                        ("vertical-align", "bottom"),
                                                        ("width", str_col_width + " !important"),
                                                        ("-webkit-width", str_col_width + " !important"),
                                                        ("-moz-width", str_col_width + " !important"),
                                                        ("-o-width", str_col_width + " !important"),
                                                        ("-ms-width", str_col_width + " !important")]),
        dict(selector="tr td:nth-child(n+" + str(num_of_cols - N + 2) + ")", props=[("text-align", "right"),
                                                        ("width", str_col_width + " !important"),
                                                        ("-webkit-width", str_col_width + " !important"),
                                                        ("-moz-width", str_col_width + " !important"),
                                                        ("-o-width", str_col_width + " !important"),
                                                        ("-ms-width", str_col_width + " !important")])
        ]


    styler = table.round(2).style.set_uuid(table_id).apply(set_cell_colors,subset=gradient_cols, axis=1).set_table_styles(styles)
    html_string = styler.render()
    for i, r in enumerate(gradient_cols):
        html_string = html_string.replace("col_heading level0 col" + str(i + num_of_cols - N) + "\" >" + str(r) + "</th>",
                                          "col_heading level0 col" + str(i + num_of_cols - N) + "\" ><div class=\"verticalTableHeader\" style=\"width:" + str_col_width + "\" >" + str(r) + "</div></th>")

    
    script1 = '<script type="text/javascript" charset="utf8" src=' \
              '"' + join(os.path.dirname(os.path.abspath(__file__)), "style/table_0.js") + '"' \
              '></script>'
    script2 = '<script type="text/javascript" charset="utf8" src=' \
              '"' + join(os.path.dirname(os.path.abspath(__file__)), "style/table_1.js") + '"' \
              '></script>'
    script3 = '<script> $(function(){$("#T_' + table_id + '").dataTable({"iDisplayLength": 50}); })</script>'

    # write combined html code
    with open(path, 'w') as file_out:
        file_out.write(title + page_css_string + html_string + script1 + script2 + script3 + table_css_string)

