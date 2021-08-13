USE `DESIRE` ;

Create table 3DStructures (
    3D_structure_id INT NOT NULL,
    StructureName varchar(50),
    PRIMARY KEY (3D_structure_id)
);

  
Create table Interactions (
    interactions_id INT NOT NULL,
    residue_i INT,
    residue_j INT,
    bp_type varchar(5),
    bp_group varchar(13),
    3D_structure_id INT NOT NULL,
    PRIMARY KEY (interactions_id),
    FOREIGN KEY (3D_structure_id) REFERENCES 3DStructures(3D_structure_id)
);


Create table ChainList (
    ChainList_id INT NOT NULL,
    3D_structure_id INT NOT NULL,
    polymer_id INT NOT NULL,
    ChainName varchar(3),
    PRIMARY KEY (ChainList_id),
    FOREIGN KEY (3D_structure_id) REFERENCES 3DStructures(3D_structure_id),
    FOREIGN KEY (polymer_id) REFERENCES Polymer_Data(PData_id)
);

Create table TextLabels (
    TextLabel_id INT NOT NULL,
    LabelText varchar(500),
    X DOUBLE(8,3),
    Y DOUBLE(8,3),
    Font varchar(100),
    Font_Size DOUBLE(6,3),
    Fill Char(50),
    secondary_structure_id INT NOT NULL,
    PRIMARY KEY (TextLabel_id),
    FOREIGN KEY (secondary_structure_id) REFERENCES SecondaryStructures(SecStr_id)
);

Create table StructuralData2 (
    map_index INT,
    Domain_RN varchar(4),
    Domain_AN int(1),
    Domains_Color INT(1),
    Helix_Num varchar(4),
    Helix_Color INT(1),
    secondary_structure_id INT,
    FOREIGN KEY (secondary_structure_id) REFERENCES SecondaryStructures(SecStr_id)
);

Create table LineLabels (
    LineLabel_id INT NOT NULL,
    X1 DOUBLE(8,3),
    Y1 DOUBLE(8,3),
    X2 DOUBLE,
    Y2 DOUBLE(8,3),
    Fill char(7),
    Stroke char(7),
    StrokeWidth DOUBLE(8,3),
    StrokeLineJoin char(5),
    StrokeMiterLimit DOUBLE(6,3),
    secondary_structure_id INT NOT NULL,
    PRIMARY KEY(LineLabel_id),
    FOREIGN KEY (secondary_structure_id) REFERENCES SecondaryStructures(SecStr_id)
);

Create table Default3DStructure (
    secondary_structure_id INT NOT NULL,
    3D_structure_id INT NOT NULL,
    FOREIGN KEY (secondary_structure_id) REFERENCES SecondaryStructures(SecStr_id),
    FOREIGN KEY (3D_structure_id) REFERENCES 3DStructures(3D_structure_id)
);

CREATE table Secondary_Tertiary (
    secondary_tertiary_id INT NOT NULL,
    secondary_structure_id INT NOT NULL,
    3D_structure_id INT NOT NULL,
    PRIMARY KEY (secondary_tertiary_id),
    FOREIGN KEY (secondary_structure_id) REFERENCES SecondaryStructures(SecStr_id),
    FOREIGN KEY (3D_structure_id) REFERENCES 3DStructures(3D_structure_id)
);


Create table StructDataMenuDetails (
    struct_data_id INT,
    StructDataName varchar(100),
    VariableName varchar(15),
    ColorList varchar(12),
    IndexMode varchar(5),
    ExtraArg varchar(9),
    Description varchar(255),
    HelpLink varchar(45),
    PRIMARY KEY (struct_data_id)
);

Create table StructDataMenu (
    StructDataMenu_id INT,
    3D_structure_id INT NOT NULL,
    struct_data_id INT NOT NULL,
    PRIMARY KEY (StructDataMenu_id),
    FOREIGN KEY (3D_structure_id) REFERENCES 3DStructures(3D_structure_id),
    FOREIGN KEY (struct_data_id) REFERENCES StructDataMenuDetails(struct_data_id)
);

Create table StructuralData3 (
    map_index INT(4),
    Value FLOAT, 
    struct_data_id INT NOT NULL,
    3D_structure_id INT NOT NULL,
    FOREIGN KEY (struct_data_id) REFERENCES StructDataMenuDetails(struct_data_id),
    FOREIGN KEY (3D_structure_id) REFERENCES 3DStructures(3D_structure_id)
);


/*Additional table to make queries easier*/
CREATE TABLE IF NOT EXISTS `DESIRE`.`Polymer_Alignments` (
  `PData_id` INT NOT NULL,
  `Aln_id` INT NOT NULL,
  PRIMARY KEY (`PData_id`, `Aln_id`),
  INDEX `alignment_fk_idx` (`Aln_id` ASC),
  CONSTRAINT `polymer_fk`
    FOREIGN KEY (`PData_id`)
    REFERENCES `DESIRE`.`Polymer_Data` (`PData_id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `alignment_fk`
    FOREIGN KEY (`Aln_id`)
    REFERENCES `DESIRE`.`Alignment` (`Aln_id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

/*Table for PDB topology viewer ids*/
CREATE TABLE IF NOT EXISTS `DESIRE`.`TopologyViewerEntities` (
  `chain_id` INT NOT NULL,
  `entity_id` INT NOT NULL,
  INDEX `chainlist_fk_idx` (`chain_id` ASC),
  PRIMARY KEY (`chain_id`),
  CONSTRAINT `chainlist_fk`
    FOREIGN KEY (`chain_id`)
    REFERENCES `DESIRE`.`ChainList` (`ChainList_id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB


/*And fill it up with data*/
INSERT INTO Polymer_Alignments (PData_id, Aln_id)
SELECT DESIRE.Polymer_Data.PData_id, DESIRE.Alignment.Aln_id FROM DESIRE.Aln_Data
INNER JOIN DESIRE.Alignment ON DESIRE.Aln_Data.aln_id = DESIRE.Alignment.Aln_id
INNER JOIN DESIRE.Residues ON DESIRE.Aln_Data.res_id = DESIRE.Residues.resi_id
INNER JOIN DESIRE.Polymer_Data ON DESIRE.Residues.PolData_id = DESIRE.Polymer_Data.PData_id
GROUP BY DESIRE.Polymer_Data.PData_id, DESIRE.Alignment.Aln_id;
