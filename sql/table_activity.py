import sqlite3

# SQLite Database Connection
conn = sqlite3.connect("chembl.db")
cursor = conn.cursor()

def create():
    # Create 'activity' Table
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS activity (
            activity_id INTEGER PRIMARY KEY,
            activity_comment TEXT,
            assay_chembl_id TEXT,
            assay_description TEXT,
            assay_type TEXT,
            assay_variant_accession TEXT,
            assay_variant_mutation TEXT,
            bao_endpoint TEXT,
            bao_format TEXT,
            bao_label TEXT,
            canonical_smiles TEXT,
            data_validity_comment TEXT,
            data_validity_description TEXT,
            document_chembl_id TEXT,
            document_journal TEXT,
            document_year INTEGER,
            molecule_chembl_id TEXT,
            molecule_pref_name TEXT,
            parent_molecule_chembl_id TEXT,
            pchembl_value REAL,
            potential_duplicate BOOLEAN,
            qudt_units TEXT,
            record_id INTEGER,
            relation TEXT,
            src_id INTEGER,
            standard_flag INTEGER,
            standard_relation TEXT,
            standard_text_value TEXT,
            standard_type TEXT,
            standard_units TEXT,
            standard_upper_value REAL,
            standard_value REAL,
            target_chembl_id TEXT,
            target_organism TEXT,
            target_pref_name TEXT,
            target_tax_id INTEGER,
            text_value TEXT,
            toid TEXT,
            type TEXT,
            units TEXT,
            uo_units TEXT,
            upper_value REAL,
            value REAL
        )
    """)
    # Create 'activity__activity_properties' Table
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS activity__activity_properties (
            id INTEGER PRIMARY KEY AUTOINCREMENT,  -- Unique ID for each property
            activity_id INTEGER,                   -- Foreign key linking to activity table
            comments TEXT,
            relation TEXT,
            result_flag INTEGER,
            standard_relation TEXT,
            standard_text_value TEXT,
            standard_type TEXT,
            standard_units TEXT,
            standard_value REAL,
            text_value TEXT,
            type TEXT,
            units TEXT,
            value REAL,
            FOREIGN KEY (activity_id) REFERENCES activity(activity_id) ON DELETE CASCADE
        )
    """)
    # Create 'activity__action_type' Table
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS activity__action_type (
            id INTEGER PRIMARY KEY AUTOINCREMENT,  -- Unique ID for each property
            activity_id INTEGER,                   -- Foreign key linking to activity table
            action_type TEXT,
            description TEXT,
            parent_type TEXT,
            FOREIGN KEY (activity_id) REFERENCES activity(activity_id) ON DELETE CASCADE
        )
    """)
    # Create 'activity__ligand_efficiency' Table
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS activity__ligand_efficiency (
            id INTEGER PRIMARY KEY AUTOINCREMENT,  -- Unique ID for each property
            activity_id INTEGER,                   -- Foreign key linking to activity table
            bei REAL,
            le REAL,
            lle REAL,
            sei REAL,
            FOREIGN KEY (activity_id) REFERENCES activity(activity_id) ON DELETE CASCADE
        )
    """)
    
    conn.commit()

def insert(activity):
    try:
        cursor.execute("""
            INSERT OR IGNORE INTO activity (
                activity_id, 
                activity_comment, 
                assay_chembl_id, 
                assay_description, 
                assay_type, 
                assay_variant_accession, 
                assay_variant_mutation, 
                bao_endpoint, 
                bao_format, 
                bao_label, 
                canonical_smiles, 
                data_validity_comment, 
                data_validity_description, 
                document_chembl_id, 
                document_journal, 
                document_year,
                molecule_chembl_id, 
                molecule_pref_name, 
                parent_molecule_chembl_id, 
                pchembl_value, 
                potential_duplicate, 
                qudt_units, 
                record_id, 
                relation, 
                src_id, 
                standard_flag, 
                standard_relation, 
                standard_text_value, 
                standard_type, 
                standard_units, 
                standard_upper_value, 
                standard_value, 
                target_chembl_id, 
                target_organism, 
                target_pref_name, 
                target_tax_id, 
                text_value, 
                toid, 
                type, 
                units, 
                uo_units, 
                upper_value, 
                value
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, 
            (
                activity["activity_id"], 
                activity["activity_comment"], 
                activity["assay_chembl_id"], 
                activity["assay_description"], 
                activity["assay_type"], 
                activity["assay_variant_accession"],
                activity["assay_variant_mutation"], 
                activity["bao_endpoint"], 
                activity["bao_format"], 
                activity["bao_label"], 
                activity["canonical_smiles"],
                activity["data_validity_comment"],
                activity["data_validity_description"], 
                activity["document_chembl_id"], 
                activity["document_journal"],
                activity["document_year"], 
                activity["molecule_chembl_id"], 
                activity["molecule_pref_name"],
                activity["parent_molecule_chembl_id"], 
                activity["pchembl_value"], 
                activity["potential_duplicate"], 
                activity["qudt_units"],
                activity["record_id"], 
                activity["relation"], 
                activity["src_id"], 
                activity["standard_flag"], 
                activity["standard_relation"],
                activity["standard_text_value"], 
                activity["standard_type"], 
                activity["standard_units"], 
                activity["standard_upper_value"],
                activity["standard_value"], 
                activity["target_chembl_id"], 
                activity["target_organism"], 
                activity["target_pref_name"], 
                activity["target_tax_id"],
                activity["text_value"], 
                activity["toid"], 
                activity["type"], 
                activity["units"], 
                activity["uo_units"], 
                activity["upper_value"], 
                activity["value"]
            )
        )
        if activity_properties := activity["activity_properties"]:
            for prop in activity_properties:
                cursor.execute("""
                    INSERT INTO activity__activity_properties (
                        activity_id, comments, relation, result_flag, standard_relation, 
                        standard_text_value, standard_type, standard_units, standard_value, 
                        text_value, type, units, value
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """, (
                    activity["activity_id"], prop["comments"], prop["relation"], prop["result_flag"],
                    prop["standard_relation"], prop["standard_text_value"], prop["standard_type"],
                    prop["standard_units"], float(prop["standard_value"]) if prop["standard_value"] else None,
                    prop["text_value"], prop["type"], prop["units"], float(prop["value"]) if prop["value"] else None
                ))
        if action_type := activity["action_type"]:
            cursor.execute("""
                INSERT INTO activity__action_type (
                    activity_id, action_type, description, parent_type
                ) VALUES (?, ?, ?, ?)
            """, (
                activity["activity_id"], action_type["action_type"], action_type["description"], action_type["parent_type"]
            ))
        if ligand_efficiency := activity["ligand_efficiency"]:
            cursor.execute("""
                INSERT INTO activity__ligand_efficiency (
                    activity_id, bei, le, lle, sei
                ) VALUES (?, ?, ?, ?, ?)
            """, (
                activity["activity_id"], ligand_efficiency["bei"], ligand_efficiency["le"], ligand_efficiency["lle"], ligand_efficiency["sei"]
            ))
        
        conn.commit()
        print("Activity data inserted successfully!")
    except Exception as e:
        print("Error inserting activity data:", e)
        import json
        print(json.dumps(activity, indent=4))
        raise e

def fetch_all():
    conn = sqlite3.connect("chembl.db")
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM activity")
    return cursor.fetchall()
