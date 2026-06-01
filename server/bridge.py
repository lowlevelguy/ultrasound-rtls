#!/usr/bin/python3
"""
Binary MQTT decoder stage.

Subscribes to `rtls/position/raw` (raw binary timepos structs from ESP32),
unpacks them, and republishes as JSON to `rtls/position` for Telegraf.

Also handles the request/response handshake:
  - Listens on `rtls/position/request`
  - Queries InfluxDB for the last recorded timestamp
  - Responds on `rtls/position/response` with that timestamp (uint32, little-endian)

struct timepos {
    uint32_t timestamp;
    float pos[2];   // [lat, lon]
};
"""

import struct
import json
import logging
import paho.mqtt.client as mqtt
from paho.mqtt.client import CallbackAPIVersion
from influxdb_client import InfluxDBClient

# ── Configuration ─────────────────────────────────────────────────────────────
MQTT_BROKER    = "127.0.0.1"
MQTT_PORT      = 1883
TOPIC_RAW      = "rtls/position/raw"
TOPIC_DECODED  = "rtls/position"
TOPIC_REQUEST  = "rtls/position/request"
TOPIC_RESPONSE = "rtls/position/response"

INFLUX_URL    = "http://127.0.0.1:8086"
INFLUX_TOKEN  = "TjTVz2tOTxIP46mlUF0qpGcMg6DyZpHB0CAG69hTNMv7_eh10DvBUG3TLYtLOhv0mKe2NHMH7Kj96-lHeXBnpA=="
INFLUX_ORG    = "rtls"
INFLUX_BUCKET = "positions"
# ─────────────────────────────────────────────────────────────────────────────

STRUCT_FMT  = "<Iff"
STRUCT_SIZE = struct.calcsize(STRUCT_FMT)  # 12 bytes

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)

influx = InfluxDBClient(url=INFLUX_URL, token=INFLUX_TOKEN, org=INFLUX_ORG)
query_api = influx.query_api()


def get_last_timestamp() -> int | None:
    query = f'''
        from(bucket: "{INFLUX_BUCKET}")
          |> range(start: 0)
          |> filter(fn: (r) => r._measurement == "position" and r._field == "ts")
          |> sort(columns: ["_time"], desc: true)
          |> limit(n: 1)
    '''
    tables = query_api.query(query)
    for table in tables:
        for record in table.records:
            log.info("raw record values: %s", dict(record.values))
            raw = int(record["_value"]); return raw // 1_000_000_000
    return None


def on_connect(client, userdata, flags, reason_code, properties):
    if reason_code == 0:
        log.info("Connected. Subscribing to '%s' and '%s'", TOPIC_RAW, TOPIC_REQUEST)
        client.subscribe(TOPIC_RAW)
        client.subscribe(TOPIC_REQUEST)
    else:
        log.error("Connection failed (reason code %s)", reason_code)


def on_message(client, userdata, msg):
    if msg.topic == TOPIC_REQUEST:
        handle_request(client)
    elif msg.topic == TOPIC_RAW:
        handle_raw(client, msg)


def handle_request(client):
    """Query InfluxDB for the last timestamp and publish it to TOPIC_RESPONSE."""
    last_ts = get_last_timestamp()

    if last_ts is None:
        # Bucket is empty — tell the ESP32 to send everything it has.
        # 0 is a safe sentinel: the ESP32 should treat it as "no prior data".
        last_ts = 0
        log.info("Bucket empty, responding with ts=0")
    else:
        log.info("Last timestamp in DB: %d", last_ts)

    log.info("About to pack: %d (type: %s)", last_ts, type(last_ts))
    client.publish(TOPIC_RESPONSE, struct.pack("<I", last_ts))
    # Send back as a raw little-endian uint32 to match the ESP32's native types
    client.publish(TOPIC_RESPONSE, struct.pack("<I", last_ts))


def handle_raw(client, msg):
    """Decode a binary timepos struct and republish as JSON for Telegraf."""
    if len(msg.payload) % STRUCT_SIZE != 0:
        log.warning(
            "Unexpected payload size: %d bytes (expected %d)",
            len(msg.payload), STRUCT_SIZE,
        )
        return
    
    entry_count = len(msg.payload) // STRUCT_SIZE
    for i in range(entry_count):
        single_payload = msg.payload[i * STRUCT_SIZE : (i+1) * STRUCT_SIZE]
        timestamp, x, y = struct.unpack(STRUCT_FMT, single_payload)

        payload = json.dumps({"ts": timestamp, "x": x, "y": y})
        client.publish(TOPIC_DECODED, payload)
        log.info("ts=%d  x=%.6f  y=%.6f → republished to '%s'", timestamp, x, y, TOPIC_DECODED)


def main():
    client = mqtt.Client(CallbackAPIVersion.VERSION2)
    client.on_connect = on_connect
    client.on_message = on_message

    client.connect(MQTT_BROKER, MQTT_PORT)

    try:
        log.info("Decoder running. Press Ctrl-C to stop.")
        client.loop_forever()
    except KeyboardInterrupt:
        log.info("Shutting down.")
    finally:
        client.disconnect()
        influx.close()

main()
