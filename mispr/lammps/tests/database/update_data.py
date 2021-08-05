from first_connection import Connect
from pymongo import MongoClient
from pprint import pprint
from bson.son import SON
import asyncio

if __name__ == "__main__":
    # # Connect to local database
    connection = Connect.get_connection()

    # # Access the test database
    db = connection.test

    # # Update a single document in the inventory collection
    db.inventory.update_one({"item": "paper"},
                            {"$set": {"size.uom": "cm", "status": "P"},
                             "$currentDate": {"lastModified": True}})
    # # # Run the loop... No idea why the tutorial is having me do this.
    # # # # There is no reference to asyncio anywhere else on the page.
    # loop = asyncio.get_event_loop()
    # loop.run_until_complete(do_update_one())

    # # Update multiple documents
    db.inventory.update_many(
        {"qty": {"$lt": 50}},
        {"$set": {"size.uom": "in", "status": "P"},
         "$currentDate": {"lastModified": True}})