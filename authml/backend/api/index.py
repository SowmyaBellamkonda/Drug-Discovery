from flask import Flask, jsonify

app = Flask(__name__)

@app.route("/")
def home():
    return jsonify({"message": "Backend running on Vercel"})

def handler(request):
    return app(request.environ, lambda *args: None)
