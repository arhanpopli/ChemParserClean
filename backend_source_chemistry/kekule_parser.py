from utils import convert_mol_format
from pprint import pprint


def get_x_y_arrow(arrow_id, doc):
    arrow_json = [i for i in doc if i['id'] == arrow_id][0]
    x, y = arrow_json['coord2D']['x'], arrow_json['coord2D']['y']
    return x, y, arrow_json


def get_arrow_length(arrow_id, doc):
    """
    If the arrow was moved, then we need to recalculate
    start and end x coords. The start x coord were recalculated
    in recalculate_x_coordinates. The length of an arrow is:
        nodes[1][x] - nodes[0][x]
    """
    # Get x, y coordinates  where an arrow starts
    x_start, y, arrow_doc = get_x_y_arrow(arrow_id, doc)
    # get arrow length
    new_arrow_start_x = arrow_doc["ctab"]['nodes'][0]['coord2D']['x']
    new_arrow_end_x = arrow_doc["ctab"]['nodes'][1]['coord2D']['x']
    arrow_length = new_arrow_end_x - new_arrow_start_x
    # calculate coordinates where an arrow ends
    x_end = x_start + arrow_length
    return x_start, x_end, y


def get_above_below_arrow(arrow_id, doc):
    arrow_x_start, arrow_x_end, arrow_y = get_arrow_length(arrow_id, doc)
    items_above_arrow = []
    items_below_arrow = []

    # find items which coordinates are between x_start and x_end.
    for item in doc:

        item_x = item['coord2D']['x']
        item_y = item['coord2D']['y']

        if arrow_x_start < item_x < arrow_x_end:
            # Those items should be either above or below an arrow
            if item_y > arrow_y:
                items_above_arrow.append(item['id'])
            else:
                items_below_arrow.append(item['id'])
    return items_above_arrow, items_below_arrow


def arrow_to_latex(arrow_doc):
    if arrow_doc['reactionType'] == "normal":
        return "\\arrow{->}"
    if arrow_doc["reactionType"] == "reversible":
        return "\\arrow{<=>}"


def compose_arrow(item, ids_latex):
    # Return latex format of an arrow with above and below items

    # get id of the arrow
    key = item.keys()[0]

    # get latex of arrow
    arrow = ids_latex[key]

    # check if there are above and below items
    if item[key]['above'] and item[key]["below"]:

        above_latex = ", ".join(
            [ids_latex[el_id] for el_id in item[key]['above']]
        )

        below_latex = ", ".join(
            [ids_latex[el_id] for el_id in item[key]['below']]
        )

        if '->' in arrow:
            # arrow = "\\arrow{->[%s][%s]}" % above_latex, % below_latex
            arrow = "\\arrow{->[" + above_latex + "][" + below_latex + "]}"

        if '<=>' in arrow:
            # arrow = "\\arrow{<=>[][]}" % above_latex, % below_latex
            arrow = "\\arrow{<=>[" + above_latex + "][" + below_latex + "]}"

    elif item[key]['above'] and not item[key]["below"]:

        above_latex = ", ".join(
            [ids_latex[el_id] for el_id in item[key]['above']]
        )

        if '->' in arrow:
            arrow = "\\arrow{->[" + above_latex + "]}"

        if '<=>' in arrow:
            # arrow = "\\arrow{<=>[{}]}".format(above_latex)
            arrow = "\\arrow{<=>[" + above_latex + "]}"

    elif not item[key]['above'] and item[key]["below"]:

        below_latex = ", ".join(
            [ids_latex[el_id] for el_id in item[key]['below']]
        )

        if '->' in arrow:
            arrow = "\\arrow{->[][" + below_latex + "{}]}"

        if '<=>' in arrow:
            # arrow = "\\arrow{<=>[][{}]}".format(below_latex)
            arrow = "\\arrow{<=>[][" + below_latex + "{}]}"
    # arrow = arrow + "\n"
    return arrow


def get_symbol(section):
    # pprint(section['obj'])
    # print('OBJECT')
    if section["obj"]["__type__"] == "Kekule.Atom":
        symbol = section["obj"]["isotopeId"]
    else:
        # pprint(section['obj'])
        symbol = section["obj"]["symbol"]

    # if section["obj"]["__type__"] == "Kekule.PseudoAtom":
        # symbol = section["obj"]["symbol"]
    return symbol


def parse_children(section):
    # kekule json does not embed "()" explicitly
    # Note: "_" is prepended to a number to make it a subscript
    # for the latex format
    subformula = "("
    for child in section['obj']['sections']:
        symbol = get_symbol(child)
        subformula += symbol
        count = child['count']
        if count > 1:
            subformula += "$_" + str(count) + "$"
    subformula += ")"
    subformula += "$_" + str(section['count']) + "$"
    return subformula


def has_children(section):
    return "sections" in section['obj']


def parse_formula(formula_json):
    final_formula = ""
    for section in formula_json['sections']:
        # if an element has it's own sections, it means that
        # it has other elements in parenthesis, like Pd(PPh3)
        # where PPh3 are child elements that are parsed in
        # parse_children function
        if has_children(section):
            children = parse_children(section)
            final_formula += children
        else:
            symbol = get_symbol(section)
            count = section['count']
            final_formula += symbol
            if count > 1:
                # Note: "_" is prepended to a number to make it a subscript
                # for the latex format
                final_formula += "$_" + str(count) + "$"
    return final_formula


def molfile_to_latex(molfile, args):
    chemfig, pdflink, error = convert_mol_format(molfile, args=args)
    if error:
        return chemfig, error
    return chemfig


def text_to_latex(doc):
    text = doc['text'].strip()
    if text == "+":
        # return "\\+"
        return "\\+{1em, 1em, 2em}"

    if "%" in text:
        return text.replace("%", "\\%")

    # convert degrees of celcius to latex
    text_vec = list(text)
    for idx, letter in enumerate(text_vec):
        if letter.isdigit() and text_vec[idx + 1] == "C":
            text_vec[idx + 1] = "^\\circ$C"
    text = "".join(text_vec)
    return text


def compose_reaction(ids_latex, list_of_items):
    reaction_latex = []
    for item in list_of_items:
        # arrow with above and below items is a dict
        if isinstance(item, dict):
            arrow = compose_arrow(item, ids_latex)
            reaction_latex.append(arrow)
        else:
            reaction_latex.append(ids_latex[item])
    return "\n".join(reaction_latex)


def recalculate_x_coordinates(items):
    """
    If a user draws a reaction and then decides to move/swap any of reaction
    elements (!!! EXCEPT ARROW), then a recalculation of coordinates has to be done.
    !!! coord2D stay the same, ctab['nodes'][0]['coords2D'] are recalculated
    So the final X coordinate is calculative this way:
        ['coord2D']['x'] + ctab[nodes][0][coord2D][x]
    """
    recalculated_items = []

    for item in items:
        if "ctab" in item.keys():
            x_2 = item['ctab']['nodes'][0]['coord2D']['x']
            item['coord2D']['x'] += x_2
            recalculated_items.append(item)
        else:
            recalculated_items.append(item)
    return recalculated_items


def generate_latex(json_doc, mol_files, args):
    # get items of the reaction
    items = json_doc['root']['children']['items']
    items = recalculate_x_coordinates(items)
    # sort items based on x coordinates
    sorted_items = sorted(items,  key=lambda x: x['coord2D']['x'])
    sorted_ids = [item['id'] for item in sorted_items]

    # a list of dictionaries: key - item's id, value - item in latex format
    ids_latex = {}
    list_of_items = []
    # When parsing the arrow, above and below elements are also parsed
    # So, not to parse them again (not as a part of arrow), I need to
    # exclude them (to skip)
    items_to_skip = []

    for item in sorted_items:

        if item["__type__"] == "Kekule.Molecule":
            # Check if it's a formula, like Pd(PPh4)4
            if "formula" in item:
                ids_latex[item['id']] = parse_formula(item['formula'])
                list_of_items.append(item['id'])
            else:
                ids_latex[item['id']] =\
                    molfile_to_latex(mol_files[item['id']], args)
                list_of_items.append(item['id'])

        if item["__type__"] == "Kekule.TextBlock":
            ids_latex[item['id']] = text_to_latex(item)
            list_of_items.append(item['id'])

        if item["__type__"] == "Kekule.Glyph.HeatSymbol":
            ids_latex[item['id']] = "$\\Delta$"
            list_of_items.append(item['id'])

        if item["__type__"] == "Kekule.Glyph.AddSymbol":
            # ids_latex[item['id']] = "\\+"
            ids_latex[item['id']] = "\\+{1em, 1em, 2em}"
            list_of_items.append(item['id'])

        if item["__type__"] == "Kekule.Glyph.ReactionArrow":
            ids_latex[item['id']] = arrow_to_latex(item)
            above, below = get_above_below_arrow(item['id'], sorted_items)
            arrow_dict = {}
            arrow_dict[item['id']] = {'above': above, 'below': below}
            list_of_items.append(arrow_dict)
            concat_el = above + below
            items_to_skip += concat_el

    # remove duplicated items
    list_of_items = [i for i in list_of_items if i not in items_to_skip]
    reaction_chemfig = compose_reaction(ids_latex, list_of_items)

    # complete_latex = latex_start + scheme_start + reaction_chemfig + latex_ends
    complete_latex = latex_start + scheme_start + reaction_chemfig +\
        scheme_stop + latex_ends

    txt_latex = scheme_start + reaction_chemfig + scheme_stop
    return complete_latex, txt_latex


def update_reaction_chemfig(reaction_txt):
    return latex_start + reaction_txt + latex_ends


"""
Kekule.Molecule
Kekule.TextBlock
Kekule.Glyph.HeatSymbol
Kekule.Glyph.AddSymbol
Kekule.Glyph.ReactionArrow
"""


latex_start = r"""
\documentclass{minimal}
\usepackage{xcolor, chemfig, mol2chemfig}
\usepackage[paperheight=10cm, paperwidth=28cm]{geometry}
\usepackage[helvet]{sfmath}

\setbondoffset{1pt}
\setatomsep{2em}
\setdoublesep{3.5pt}
\setbondstyle{line width=1pt}

\begin{document}
\vspace*{\fill}
\vspace{-4pt}
\begin{center}""" + "\n"

scheme_start = "\schemestart[0, 2.0, thick]" + "\n"
# stop = "\n" + "\schemestop"
plus_symbol = "\+{1em, 1em, -2em}"
scheme_stop = "\n" + "\schemestop"
latex_ends = "\n" + "\end{center}" + "\n" \
        + "\\vspace*{\\fill}" + "\n" + "\end{document}"
# scheme_stop = "\n" + "\schemestop"
# latex_ends = "\n" + "\schemestop" + "\n" + "\end{center}" + "\n" \

        # + "\\vspace*{\\fill}" + "\n" + "\end{document}"
