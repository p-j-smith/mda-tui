from MDAnalysis import Universe
from MDAnalysis.transformations import (
    NoJump as nojump,
    TransformationBase,
    center_in_box,
    translate,
    unwrap,
    wrap,
)
from textual import on
from textual.app import ComposeResult
from textual.containers import Horizontal, Vertical
from textual.widgets import (
    ContentSwitcher,
    Input,
    Label,
    Select,
    Static,
    Switch,
)

from mda_tui.validators import AtomSelectionValidator, NumberOrNoneValidator


class WidgetWithLabel(Static):
    """Add a label to a widget"""

    def __init__(
        self,
        label: str,
        widget: Static,
        id: str | None = None,
        tooltip: str = "",
    ) -> None:
        super().__init__(id=id)
        self.label_text = label
        self.labelled_widget = widget
        self.tooltip = tooltip

    def compose(self) -> ComposeResult:
        yield Label(self.label_text)
        yield self.labelled_widget


class TransformationSelector(Vertical):
    """Widget for selecting a transformation to apply"""

    def compose(self) -> ComposeResult:
        """Create the layout for the transformation selector"""

        translate = Translate(id=Translate.id)
        center_in_box = CenterInBox(id=CenterInBox.id)
        wrap = Wrap(id=Wrap.id)
        unwrap = Unwrap(id=Unwrap.id)
        nojump = NoJump(id=NoJump.id)

        options = [
            (translate.description, translate),
            (center_in_box.description, center_in_box),
            (wrap.description, wrap),
            (unwrap.description, unwrap),
            (nojump.description, nojump),
        ]
        select: Select[TransformationBase] = Select(
            options=options,
            prompt="select transformation",
            id="transformation",
        )
        select.tooltip = "select a transformation to apply"
        yield select
        with ContentSwitcher(initial=None, id="selector"):
            yield translate
            yield center_in_box
            yield wrap
            yield unwrap
            yield nojump

    @on(Select.Changed, "#transformation")
    def select_transformation(self, event) -> None:
        selected_id = None if event.value is None else event.value.id
        self.query_one(ContentSwitcher).current = selected_id


class Translate(Vertical):
    """Widgets for setting parameters for translate transformation"""

    description = "Translate coordinates by a given vector"
    transformation = translate
    id = str(translate).removeprefix("<class '").removesuffix("'>")  # noqa: A003

    def compose(self) -> ComposeResult:
        """Create layout of parameter widgets"""

        vector = Horizontal(
            Input(
                placeholder="x",
                validators=NumberOrNoneValidator(
                    failure_description="'x' must be a number or empty",
                ),
                id="translate_x",
            ),
            Input(
                placeholder="y",
                validators=NumberOrNoneValidator(
                    failure_description="'y' must be a number or empty",
                ),
                id="translate_y",
            ),
            Input(
                placeholder="z",
                validators=NumberOrNoneValidator(
                    failure_description="'z' must be a number or empty",
                ),
                id="translate_z",
            ),
        )

        # Define tooltips
        vector_tooltip = "vector by which the all atoms will be translated"

        yield WidgetWithLabel(
            label="vector",
            widget=vector,
            id="translate_vector",
            tooltip=vector_tooltip,
        )

    @property
    def vector(self):
        x = self.query_one("#translate_x", Input).value
        y = self.query_one("#translate_y", Input).value
        z = self.query_one("#translate_z", Input).value
        x = float(x) if x else 0
        y = float(y) if y else 0
        z = float(z) if z else 0
        return [x, y, z]

    def setup_transformation(self, universe: Universe):  # noqa: ARG002
        """Initialise the transformation for a given universe"""
        return self.transformation(vector=self.vector)

    def validate(self):
        return [widget.validate(widget.value) for widget in self.query(Input)]


class CenterInBox(Vertical):
    """Widgets for setting parameters for center_in_box transformation"""

    description = "Center atoms / molecules"
    transformation = center_in_box
    id = str(center_in_box).removeprefix("<class '").removesuffix("'>")  # noqa: A003

    def compose(self) -> ComposeResult:
        """Create layout of parameter widgets"""

        ag = Input(
            placeholder="atom selection",
            validators=AtomSelectionValidator(
                failure_description="invalid atom selection",
            ),
            id="ag",
        )

        options = [
            ("centor of geometry", "geometry"),
            ("center of mass", "mass"),
        ]
        method: Select[str] = Select(
            options=options,
            prompt="select centering method",
            value="geometry",
            id="centering_method",
        )

        point = Horizontal(
            Input(
                placeholder="x",
                validators=NumberOrNoneValidator(
                    failure_description="'x' must be a number or empty",
                ),
                id="center_x",
            ),
            Input(
                placeholder="y",
                validators=NumberOrNoneValidator(
                    failure_description="'y' must be a number or empty",
                ),
                id="center_y",
            ),
            Input(
                placeholder="z",
                validators=NumberOrNoneValidator(
                    failure_description="'z' must be a number or empty",
                ),
                id="center_z",
            ),
        )
        wrap = Switch()

        # Define tooltips
        ag_tooltip = "selection string for atom group to be centered on the unit cell. If empty, will center all atoms."
        method_tooltip = "choose the method of centering on the selected AtomGroup."
        point_tooltip = "if provided, overrides the unit cell center - the coordinates of each timestep are translated so that the center of mass/geometry of the selected AtomGroup is aligned to this position instead. Leave empty to translate to the center of the unit cell."
        wrap_tooltip = "if enabled, the atoms from the selected AtomGroup will be moved to the primary unit cell before calculating the center of mass or geometry."

        yield WidgetWithLabel(label="ag", widget=ag, id="center_ag", tooltip=ag_tooltip)
        yield WidgetWithLabel(
            label="center",
            widget=method,
            id="center_method",
            tooltip=method_tooltip,
        )
        yield WidgetWithLabel(label="point", widget=point, id="center_on", tooltip=point_tooltip)
        yield WidgetWithLabel(label="wrap", widget=wrap, id="center_wrap", tooltip=wrap_tooltip)

    @property
    def selection(self):
        sel = self.query_one("#ag", Input).value
        sel = sel if sel else "all"
        return sel

    @property
    def method(self):
        return self.query_one(Select).value

    @property
    def point(self):
        x = self.query_one("#center_x", Input).value
        y = self.query_one("#center_y", Input).value
        z = self.query_one("#center_z", Input).value
        x = float(x) if x else None
        y = float(y) if y else None
        z = float(z) if z else None
        return None if None in [x, y, z] else [x, y, z]

    @property
    def wrap(self):
        return self.query_one(Switch).value

    def setup_transformation(self, universe: Universe):
        """Initialise the transformation for a given universe"""
        ag = universe.select_atoms(self.selection)
        return self.transformation(ag=ag, center=self.method, point=self.point, wrap=self.wrap)

    def validate(self):
        return [widget.validate(widget.value) for widget in self.query(Input)]


class Wrap(Vertical):
    """Widgets for setting parameters for the wrap transformation"""

    description = "Wrap atoms into the unit cell"
    transformation = wrap
    id = str(wrap).removeprefix("<class '").removesuffix("'>")  # noqa: A003

    """
    ag: Atomgroup
        Atomgroup to be wrapped in the unit cell
    compound : {'atoms', 'group', 'residues', 'segments', 'fragments'}, optional
        The group which will be kept together through the shifting process.
    """

    def compose(self) -> ComposeResult:
        """Create layout of parameter widgets"""

        ag = Input(
            placeholder="atom selection",
            validators=AtomSelectionValidator(
                failure_description="invalid atom selection",
            ),
            id="ag",
        )

        options = [
            ("individual atoms", "atoms"),
            ("selected atom group", "group"),
            ("residues", "residues"),
            ("segments", "segments"),
            ("fragments", "fragments"),
        ]
        compound: Select[str] = Select(
            options=options,
            prompt="select compound to wrap",
            value="atoms",
            id="wrap_compound",
        )

        # define tooltips
        ag_tooltip = "selection string for atom group to be wrapped into the unit cell. If empty, all atoms will be wrapped."
        compound_tooltip = (
            "select the group that will be kept together through the wrapping transformation"
        )

        yield WidgetWithLabel(label="ag", widget=ag, id="wrap_ag", tooltip=ag_tooltip)
        yield WidgetWithLabel(
            label="compound",
            widget=compound,
            id="wrap_compound",
            tooltip=compound_tooltip,
        )

    @property
    def selection(self):
        sel = self.query_one("#ag", Input).value
        sel = sel if sel else "all"
        return sel

    @property
    def compound(self):
        return self.query_one(Select).value

    def setup_transformation(self, universe: Universe):
        """Initialise the transformation for a given universe"""
        ag = universe.select_atoms(self.selection)
        return self.transformation(ag=ag, compound=self.compound)

    def validate(self):
        return [widget.validate(widget.value) for widget in self.query(Input)]


class Unwrap(Vertical):
    """Widgets for setting parameters for the unwrap transformation"""

    description = "Unwrap an atom group"
    transformation = unwrap
    id = str(unwrap).removeprefix("<class '").removesuffix("'>")  # noqa: A003

    def compose(self) -> ComposeResult:
        """Create layout of parameter widgets"""

        ag = Input(
            placeholder="atom selection",
            validators=AtomSelectionValidator(
                failure_description="invalid atom selection",
            ),
            id="ag",
        )

        # define tooltips
        ag_tooltip = "selection string for atom group to be unwrapped into the unit cell. If empty, all fragments will be unwrapped."

        yield WidgetWithLabel(label="ag", widget=ag, id="unwrap_ag", tooltip=ag_tooltip)

    @property
    def selection(self):
        sel = self.query_one("#ag", Input).value
        sel = sel if sel else "all"
        return sel

    def setup_transformation(self, universe: Universe):
        """Initialise the transformation for a given universe"""
        ag = universe.select_atoms(self.selection)
        return self.transformation(ag=ag)

    def validate(self):
        return [widget.validate(widget.value) for widget in self.query(Input)]


class NoJump(Vertical):
    """Widgets for setting parameters for nojump transformation"""

    description = "Prevent jumps across periodic boundaries"
    transformation = nojump
    id = str(nojump).removeprefix("<class '").removesuffix("'>")  # noqa: A003

    def compose(self) -> ComposeResult:
        """Create layout of parameter widgets"""

        check_continuity = Switch(value=True)

        # Define tooltips
        check_continuity_tooltip = "check the trajectory step size is not more than 1"

        yield WidgetWithLabel(
            label="check continuity",
            widget=check_continuity,
            id="check_continuity",
            tooltip=check_continuity_tooltip,
        )

    @property
    def check_continuity(self):
        return self.query_one(Switch).value

    def setup_transformation(self, universe: Universe):  # noqa: ARG002
        """Initialise the transformation for a given universe"""
        return self.transformation(check_continuity=self.check_continuity)

    def validate(self):
        return [widget.validate(widget.value) for widget in self.query(Input)]
