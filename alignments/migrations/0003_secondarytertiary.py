# Generated by Django 2.2.5 on 2020-03-30 14:34

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('alignments', '0002_chainlist_threedstructures'),
    ]

    operations = [
        migrations.CreateModel(
            name='SecondaryTertiary',
            fields=[
                ('secondary_tertiary_id', models.IntegerField(primary_key=True, serialize=False)),
            ],
            options={
                'db_table': 'Secondary_Tertiary',
                'managed': False,
            },
        ),
    ]