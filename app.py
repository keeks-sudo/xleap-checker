"""
XLEAP Index Color Balance Checker

A web-based tool to validate index combinations for
Illumina XLEAP-SBS chemistry (NovaSeq X Plus, NextSeq 1000/2000).
"""

from flask import Flask, render_template, request, flash, redirect, url_for
from core.validator import (
    validate_from_text,
    parse_indexes,
    get_color_balance_table,
    XLEAP_COLORS,
)

app = Flask(__name__)
app.secret_key = 'dev-secret-key-change-in-production'


@app.route('/')
def index():
    """Home page with index input form."""
    return render_template('index.html')


@app.route('/check', methods=['POST'])
def check():
    """Process indexes and show results."""
    index_text = request.form.get('indexes', '').strip()

    if not index_text:
        flash('Please enter at least one index sequence.', 'error')
        return redirect(url_for('index'))

    # Parse and validate
    indexes = parse_indexes(index_text)

    if not indexes:
        flash('No valid index sequences found. Use only A, C, G, T characters.', 'error')
        return redirect(url_for('index'))

    # Run validation
    result = validate_from_text(index_text)
    table = get_color_balance_table(indexes)

    return render_template(
        'results.html',
        result=result,
        table=table,
        indexes=indexes,
        colors=XLEAP_COLORS,
    )


@app.route('/about')
def about():
    """About page explaining the tool."""
    return render_template('about.html')


@app.route('/help')
def help_page():
    """Help page with XLEAP chemistry explanation."""
    return render_template('help.html')


if __name__ == '__main__':
    app.run(debug=True, port=5001)
