import redis
import json

r = redis.Redis(host='localhost', port=6379, decode_responses=True)

def cache_get(key):
    value = r.get(key)
    if value:
        return json.loads(value)
    return None

def cache_set(key, data, ex=3600):
    r.set(key, json.dumps(data), ex=ex)
