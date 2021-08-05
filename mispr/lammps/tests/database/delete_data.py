from first_connection import Connect
from pymongo import MongoClient
from pprint import pprint
from bson.son import SON

if __name__ == "__main__":
    # # Connect to local database
    connection = Connect.get_connection()

    # # Access the test database
    db = connection.test

    # # Delete single document
    db.inventory.delete_one({"status": "D"})

    # # Delete multiple documents
    db.inventory.delete_many({"status": "A"})