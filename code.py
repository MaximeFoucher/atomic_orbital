import numpy as np
from scipy.special import genlaguerre, lpmv
from math import factorial
import plotly.graph_objects as go
from dash import Dash, dcc, html, Input, Output, callback

a_0 = 1

def radial(n, l, r):
    rho = 2 * r / (n * a_0)
    norm = np.sqrt((2 / (n * a_0))**3 * factorial(n-l-1) / (2 * n * factorial(n+l)))
    return norm * np.exp(-rho / 2) * rho**l * genlaguerre(n-l-1, 2*l+1)(rho)

def angular(l, m, theta, phi):
    norm = np.sqrt((2*l+1) / (4*np.pi) * factorial(l-abs(m)) / factorial(l+abs(m)))
    P = lpmv(abs(m), l, np.cos(theta))
    if m == 0:
        return norm * P
    elif m > 0:
        return norm * P * np.sqrt(2) * np.cos(m * phi)
    else:
        return norm * P * np.sqrt(2) * np.sin(abs(m) * phi)

def psi_squared(x, y, z, n, l, m):
    r     = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(np.clip(z / (r + 1e-10), -1, 1))
    phi   = np.arctan2(y, x)
    return (radial(n, l, r) * angular(l, m, theta, phi))**2

def generate_cloud(n, l, m, N=1000000):
    scale = n**2 * 3
    x = np.random.uniform(-scale, scale, N)
    y = np.random.uniform(-scale, scale, N)
    z = np.random.uniform(-scale, scale, N)
    prob = psi_squared(x, y, z, n, l, m)
    prob_norm = prob / prob.max()
    return x, y, z, prob_norm

def make_figure(n, l, m):
    x, y, z, prob = generate_cloud(n, l, m)
    mask = prob > 0.001
    labels = ['s', 'p', 'd', 'f']
    lim = n**2 * 3

    fig = go.Figure(data=go.Scatter3d(
        x=x[mask], y=y[mask], z=z[mask],
        mode='markers',
        marker=dict(
            size=1.5,
            color=prob[mask],
            colorscale='Hot',
            opacity=0.6,
            colorbar=dict(title='|ψ|²', thickness=15),
            cmin=0, cmax=1,
        )
    ))
    fig.update_layout(
        title=dict(
            text=f"Orbitale {n}{labels[l]} | m={m} | E = {-13.6/n**2:.3f} eV",
            font=dict(size=16, color='white')
        ),
        scene=dict(
            xaxis=dict(title='x (a₀)', range=[-lim, lim]),
            yaxis=dict(title='y (a₀)', range=[-lim, lim]),
            zaxis=dict(title='z (a₀)', range=[-lim, lim]),
            bgcolor='rgb(10, 10, 20)',
            xaxis_backgroundcolor='rgb(15, 15, 30)',
            yaxis_backgroundcolor='rgb(15, 15, 30)',
            zaxis_backgroundcolor='rgb(15, 15, 30)',
        ),
        paper_bgcolor='rgb(15, 15, 30)',
        font_color='white',
        margin=dict(l=0, r=0, t=50, b=0),
        height=650,
        uirevision='constant',   # conserve la rotation entre updates
    )
    return fig

# --- App Dash ---
app = Dash(__name__)

app.layout = html.Div(
    style={'backgroundColor': '#0a0a14', 'padding': '20px', 'fontFamily': 'monospace'},
    children=[
        html.H2("Orbitales atomiques — Hydrogène",
                style={'color': 'white', 'textAlign': 'center'}),

        # Sliders
        html.Div(style={'display': 'flex', 'gap': '40px', 'justifyContent': 'center',
                        'marginBottom': '20px'},
            children=[
                html.Div([
                    html.Label("n  (niveau principal)", style={'color': '#aaa', 'fontSize': '13px'}),
                    dcc.Slider(id='slider-n', min=1, max=5, step=1, value=3,
                               marks={i: str(i) for i in range(1, 6)},
                               tooltip={"placement": "bottom"})
                ], style={'width': '250px'}),

                html.Div([
                    html.Label("l  (sous-couche)", style={'color': '#aaa', 'fontSize': '13px'}),
                    dcc.Slider(id='slider-l', min=0, max=4, step=1, value=2,
                               marks={i: str(i) for i in range(5)},
                               tooltip={"placement": "bottom"})
                ], style={'width': '250px'}),

                html.Div([
                    html.Label("m  (orientation)", style={'color': '#aaa', 'fontSize': '13px'}),
                    dcc.Slider(id='slider-m', min=-4, max=4, step=1, value=0,
                               marks={i: str(i) for i in range(-4, 5)},
                               tooltip={"placement": "bottom"})
                ], style={'width': '300px'}),
            ]
        ),

        # Info quantique
        html.Div(id='info-box', style={
            'color': '#7cd4ff', 'textAlign': 'center',
            'fontSize': '13px', 'marginBottom': '10px'
        }),

        # Graphe
        dcc.Graph(id='orbital-graph', config={'scrollZoom': True}),
    ]
)

@callback(
    Output('orbital-graph', 'figure'),
    Output('info-box', 'children'),
    Output('slider-l', 'max'),       # adapter le max de l selon n
    Output('slider-m', 'min'),       # adapter min/max de m selon l
    Output('slider-m', 'max'),
    Input('slider-n', 'value'),
    Input('slider-l', 'value'),
    Input('slider-m', 'value'),
)
def update(n, l, m):
    # Contraintes quantiques
    l = min(l, n - 1)
    m = int(np.clip(m, -l, l))

    labels = ['s', 'p', 'd', 'f', 'g']
    info = (f"Règles : l ∈ [0, n−1] = [0, {n-1}]  |  "
            f"m ∈ [−l, l] = [{-l}, {l}]  |  "
            f"Sous-couche : {n}{labels[l]}  |  "
            f"Dégénérescence : {2*l+1} états")

    fig = make_figure(n, l, m)
    return fig, info, n-1, -l, l

if __name__ == '__main__':
    app.run(debug=False)